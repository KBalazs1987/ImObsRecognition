setwd("D:/CloudStation/mygit/ImObsRecognition/")
library(protr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(Rfast)

proteome = readFASTA("D:/CloudStation/mygit_EXT/ImObsRecognition/objects/uniprot-proteome_UP000005640+reviewed_yes.fasta")
aafreq = as.data.frame(Table(unlist(strsplit(unlist(proteome, use.names = F), ""), use.names = F)) / length(unlist(strsplit(unlist(proteome, use.names = F), ""), use.names = F)))
colnames(aafreq) = "freq_in_human_proteome"
rm(proteome)
aafreq["U",]

load("ligands_nonhuman_d1d2")
aacomp_nim = ligands_nonhuman %>% filter(immunogenicity == 0) %>% pull(tcem) %>% strsplit(split = "") %>% unlist(use.names = F) %>% factor(levels = rownames(AABLOSUM62)) %>% table() %>% as.numeric()
aacomp_im = ligands_nonhuman %>% filter(immunogenicity == 1) %>% pull(tcem) %>% strsplit(split = "") %>% unlist(use.names = F) %>% factor(levels = rownames(AABLOSUM62)) %>% table() %>% as.numeric()
aacomp = data.frame(aa = rownames(AABLOSUM62), nim = aacomp_nim, im = aacomp_im)
ftres = t(sapply(rownames(AABLOSUM62), function(y) {
  ft = fisher.test(matrix(data = c(aacomp$im[aacomp$aa == y],aacomp$nim[aacomp$aa == y],
                                   sum(aacomp$im[aacomp$aa != y]),sum(aacomp$nim[aacomp$aa != y])), nrow = 2))
  c(ft$estimate, ft$p.value)
}))
colnames(ftres) = c("or", "p")
rm(aacomp, aacomp_im, aacomp_nim, ligands_nonhuman)


#Excluding amino acids
resdf_exclaa = data.frame(aa = rownames(AABLOSUM62), 
                          freqpenta_med_im = NA, freqpenta_med_nim = NA, freqpenta_wilp = NA,
                          expr_fishor1 = NA, expr_fishp1 = NA, expr_fishor2 = NA, expr_fishp2 = NA,
                          thymopscore_fishor = NA, thymopscore_fishp = NA,
                          immunopscore_fishor = NA, immunopscore_fishp = NA)

for(i in 1:nrow(resdf_exclaa)) {
  load("ligands_nonhuman_d1d2")
  ligands_nonhuman %<>% filter(!grepl(resdf_exclaa$aa[i], tcem))
  resdf_exclaa$freqpenta_med_im[i] = ligands_nonhuman %>% filter(immunogenicity == 1) %>% pull(freqpenta_tcem) %>% median()
  resdf_exclaa$freqpenta_med_nim[i] = ligands_nonhuman %>% filter(immunogenicity == 0) %>% pull(freqpenta_tcem) %>% median()
  resdf_exclaa$freqpenta_wilp[i] = wilcox.test(freqpenta_tcem ~ immunogenicity, ligands_nonhuman)$p.value
  #EXPRESSION
  ligands_nonhuman$pentamer_expression_group = cut(x = ligands_nonhuman$pentamer_expression, breaks = c(min(ligands_nonhuman$pentamer_expression,na.rm=T),quantile(ligands_nonhuman$pentamer_expression,c(0.15,0.75),na.rm=T),max(ligands_nonhuman$pentamer_expression,na.rm=T)), labels = F, include.lowest = T, right = F)
  ft_lm = fisher.test(matrix(c(
    ligands_nonhuman %>% filter(immunogenicity == 0, pentamer_expression_group == 2) %>% nrow(),
    ligands_nonhuman %>% filter(immunogenicity == 1, pentamer_expression_group == 2) %>% nrow(), 
    ligands_nonhuman %>% filter(immunogenicity == 0, pentamer_expression_group == 1) %>% nrow(),
    ligands_nonhuman %>% filter(immunogenicity == 1, pentamer_expression_group == 1) %>% nrow()),
    nrow = 2, ncol = 2))
  resdf_exclaa$expr_fishor1[i] = ft_lm$estimate
  resdf_exclaa$expr_fishp1[i] = ft_lm$p.value
  ft_hm = fisher.test(matrix(c(
    ligands_nonhuman %>% filter(immunogenicity == 0, pentamer_expression_group == 2) %>% nrow(),
    ligands_nonhuman %>% filter(immunogenicity == 1, pentamer_expression_group == 2) %>% nrow(), 
    ligands_nonhuman %>% filter(immunogenicity == 0, pentamer_expression_group == 3) %>% nrow(),
    ligands_nonhuman %>% filter(immunogenicity == 1, pentamer_expression_group == 3) %>% nrow()),
    nrow = 2, ncol = 2))
  resdf_exclaa$expr_fishor2[i] = ft_hm$estimate
  resdf_exclaa$expr_fishp2[i] = ft_hm$p.value
  #THYMOPSCORE
  ligands_nonhuman$thymopscore_group = cut(x = ligands_nonhuman$thymopscore, breaks = quantile(ligands_nonhuman$thymopscore,c(0,0.25,1),na.rm = T), labels = F, include.lowest = T, right = F)
  ft_thymo = fisher.test(matrix(c(
    ligands_nonhuman %>% filter(immunogenicity == 0, thymopscore_group == 2) %>% nrow(),
    ligands_nonhuman %>% filter(immunogenicity == 1, thymopscore_group == 2) %>% nrow(), 
    ligands_nonhuman %>% filter(immunogenicity == 0, thymopscore_group == 1) %>% nrow(),
    ligands_nonhuman %>% filter(immunogenicity == 1, thymopscore_group == 1) %>% nrow()),
    nrow = 2, ncol = 2))
  resdf_exclaa$thymopscore_fishor[i] = ft_thymo$estimate
  resdf_exclaa$thymopscore_fishp[i] = ft_thymo$p.value
  #IMMUNOPSCORE
  ligands_nonhuman$immunopscore_group = cut(x = ligands_nonhuman$immunopscore, breaks = quantile(ligands_nonhuman$immunopscore,c(0,0.25,1),na.rm = T), labels = F, include.lowest = T, right = F)
  ft_immuno = fisher.test(matrix(c(
    ligands_nonhuman %>% filter(immunogenicity == 0, immunopscore_group == 2) %>% nrow(),
    ligands_nonhuman %>% filter(immunogenicity == 1, immunopscore_group == 2) %>% nrow(), 
    ligands_nonhuman %>% filter(immunogenicity == 0, immunopscore_group == 1) %>% nrow(),
    ligands_nonhuman %>% filter(immunogenicity == 1, immunopscore_group == 1) %>% nrow()),
    nrow = 2, ncol = 2))
  resdf_exclaa$immunopscore_fishp[i] = ft_immuno$p.value
  resdf_exclaa$immunopscore_fishor[i] = ft_immuno$estimate
  rm(ft_lm, ft_hm, ft_thymo, ft_immuno, ligands_nonhuman)
}
rm(i)

resdf = cbind.data.frame(freq_in_human_proteome = aafreq[rownames(aafreq) != "U",], ftres)
resdf = cbind.data.frame(resdf, resdf_exclaa[match(rownames(resdf), resdf_exclaa$aa),])
resdf = resdf[,c(4,1:3,5:15)]


