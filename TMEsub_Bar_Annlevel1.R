library(Hmisc)
library(ggplot2)
library(ggsci)
library(ggthemes)

matt = openxlsx::read.xlsx("TableS2.xlsx")

mat = matt[,-3:-1]
rownames(mat) = matt$DonorID
corl = rcorr(t(mat))
cor = corl$r
cor[which(corl$P > 0.05)] = 0
hc = hclust(as.dist(1-cor), method = "average")

matt$DonorID = factor(matt$DonorID, levels = hc$labels[hc$order])

stat = reshape2::melt(mat )
stat$Var1 = factor(stat$Var1, levels = hc$labels[hc$order])

levle1_col = c(Treg = '#bec1d4', `CD4+ T` = '#d6bcc0', `CD8+ T` = '#bb7784', NK = '#8dd593', 
               NKT = '#f0b98d', ILC = '#f6c4e1', B = '#023fa5', `Tn` = '#7D58B9',
               PC = '#11c638', Monocyte ='#F0E442',DC = '#ef9708',Mast = '#ead3c6', Neutrophil = '#ade87c', 
               Macrophage = '#FEC260', Endocrine = '#8e063b', EC = '#d33f6a', MDSC = '#e6afb9', 
               Schwann = '#0d6c0d' ,PSC = '#0fcfc0', Fibroblast = '#7382BC', SMC = '#e07b91', Adipocyte = '#F7C394')

p1 = ggplot(stat, aes(x = Var1, y = value, fill = Var2))+
  geom_bar(stat = "identity", position = "fill")+
  scale_fill_manual(values = level1_col)+
  theme_pander()+
  facet_grid(type~., scales = "free")+ xlab(NULL)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0), 
        axis.ticks.x = element_blank())

pann = ggplot(matt, aes(x = DonorID, y = 1, fill = TME_sub))+
  geom_tile()+ xlab(NULL)+
  theme_pander()+
  theme(legend.position = "top",
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())+
  scale_fill_aaas()

pann2 = ggplot(matt, aes(x = DonorID, y = 1, fill = type))+
  geom_tile()+ xlab(NULL)+
  theme_pander()+
  theme(legend.position = "top",
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

library(patchwork)
pann/pann2/p1
ggsave("Fig2A.pdf", height = 7, width = 12)
