import warnings
from typing import Optional

import scanpy as sc
import scvi
import numpy as np
from anndata import AnnData


def scale_scvi(
    adata: AnnData,
    out_path: Optional[str] = None,
    n_hvg: int = 5000,
    hvg_batch: bool = True,
    batch_key_hvg: str = "DonorID",
    batch_key_scvi: str = "DonorID",
    tech_key: str = "Tech",
    celltype_key: str = "Ann_Level1",
    run_cellhint: bool = True,
    resolution: float = 2.0,
    n_layers: int = 2,
    n_latent: int = 30,
    gene_likelihood: str = "nb",
    n_neighbors: int = 20,
    accelerator: Optional[str] = None,  # "gpu" / "cpu" / None(auto)
) -> AnnData:
    """
    Run an scVI workflow: QC -> HVG -> scVI latent -> neighbors -> UMAP (and optional CellHint).

    Parameters
    ----------
    adata : AnnData
        Input AnnData (counts in .X or .layers['counts']).
    out_path : str, optional
        If provided, write the processed AnnData to this .h5ad path.
    n_hvg : int
        Number of highly variable genes to select (Seurat v3 flavor).
    hvg_batch : bool
        If True, select HVGs in a batch-aware manner using `batch_key_hvg`.
    batch_key_hvg : str
        .obs column used for batch-aware HVG selection.
    batch_key_scvi : str
        .obs column used as scVI batch key.
    tech_key : str
        .obs column representing technology/platform (used by CellHint).
    celltype_key : str
        .obs column with coarse cell types (used by CellHint).
    run_cellhint : bool
        If True and cellhint is installed, run cellhint-based integration/denoising.
    resolution : float
        Clustering resolution if you later run Leiden (kept for compatibility).
    n_layers : int
        Number of hidden layers in scVI encoder/decoder.
    n_latent : int
        Dimensionality of scVI latent space.
    gene_likelihood : {"nb","zinb","poisson"}
        Likelihood model in scVI.
    n_neighbors : int
        Neighbors for graph construction (used by UMAP).
    accelerator : {"gpu","cpu",None}
        Training device; None lets scvi auto-detect.

    Returns
    -------
    AnnData
        A *copy* of the input with:
        - .layers['counts'] set
        - HVG-selected var (if subset happened)
        - .obsm['X_scVI'] latent
        - neighbors graph + UMAP embedding
        - .raw preserved as the original (preprocessing) snapshot
    """
    adata = adata.copy()  # work on a copy, keep caller's AnnData intact

    # --- Preserve a snapshot of the original object in .raw
    adata.raw = adata.copy()

    # --- Ensure counts layer
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X

    # --- Minimal QC (keep super light; avoid over-aggressive filtering here)
    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.filter_cells(adata, min_counts=1)

    # --- HVG selection
    if hvg_batch:
        if batch_key_hvg not in adata.obs:
            warnings.warn(
                f"batch_key_hvg='{batch_key_hvg}' not found in .obs; "
                "falling back to non-batch HVG selection."
            )
            hvg_batch = False

    if hvg_batch:
        sc.pp.highly_variable_genes(
            adata,
            flavor="seurat_v3",
            n_top_genes=n_hvg,
            layer="counts",
            batch_key=batch_key_hvg,
            subset=True,
        )
    else:
        sc.pp.highly_variable_genes(
            adata, flavor="seurat_v3", n_top_genes=n_hvg, layer="counts", subset=True
        )

    # --- scVI setup and training
    if batch_key_scvi not in adata.obs:
        warnings.warn(
            f"batch_key_scvi='{batch_key_scvi}' not found in .obs; using no batch key."
        )
        scvi.model.SCVI.setup_anndata(adata, layer="counts")
    else:
        scvi.model.SCVI.setup_anndata(
            adata, layer="counts", batch_key=batch_key_scvi
        )

    model = scvi.model.SCVI(
        adata,
        n_layers=n_layers,
        n_latent=n_latent,
        gene_likelihood=gene_likelihood,
    )

    # Auto-choose accelerator if not specified
    if accelerator is None:
        try:
            import torch

            accelerator = "gpu" if torch.cuda.is_available() else "cpu"
        except Exception:
            accelerator = "cpu"

    model.train(accelerator=accelerator)

    # --- Store latent representation
    SCVI_LATENT_KEY = "X_scVI"
    adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

    # --- Optional: CellHint integration
    if run_cellhint:
        try:
            import cellhint  # type: ignore

            if tech_key in adata.obs and celltype_key in adata.obs:
                cellhint.integrate(
                    adata,
                    batch=tech_key,
                    cell_type=celltype_key,
                    use_rep=SCVI_LATENT_KEY,
                )
            else:
                warnings.warn(
                    f"cellhint skipped: '{tech_key}' or '{celltype_key}' not in .obs."
                )
        except Exception as e:
            warnings.warn(f"cellhint not run or failed gracefully: {e}")

    # --- Graph + UMAP (use latent as representation)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=SCVI_LATENT_KEY)
    sc.tl.umap(adata)

    # NOTE: do NOT remove neighbors from .uns; UMAP depends on it
    # If you really need to drop it later for size reasons, do it after writing.

    # --- Optionally save
    if out_path is not None:
        if not out_path.endswith(".h5ad"):
            warnings.warn("out_path does not end with .h5ad; appending suffix.")
            out_path = f"{out_path}.h5ad"
        adata.write_h5ad(out_path)

    return adata
