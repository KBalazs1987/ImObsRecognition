setwd("D:/CloudStation/mygit/ImObsRecognition/")
library(fastmatch)
library(pbapply)
library(parallel)

aLign = function(seq_1, seq_2, mtx = protr::AABLOSUM62, aas = colnames(protr::AABLOSUM62)) {
  sum(mtx[cbind(fmatch(unlist(strsplit(seq_1, "")), aas), fmatch(unlist(strsplit(seq_2, "")), aas))])
}

load("ligands_nonhuman_d1d2")

ligands_nonhuman$group = 0
ligands_nonhuman$group[ligands_nonhuman$freqpenta_tcem < 4 & ligands_nonhuman$pentamer_expression < quantile(ligands_nonhuman$pentamer_expression,0.15,na.rm = T) & ligands_nonhuman$thymopscore < quantile(ligands_nonhuman$thymopscore,0.25,na.rm = T)] = 1
ligands_nonhuman$group[ligands_nonhuman$freqpenta_tcem >= 4 & (ligands_nonhuman$pentamer_expression >= quantile(ligands_nonhuman$pentamer_expression,0.15,na.rm = T) & ligands_nonhuman$pentamer_expression <= quantile(ligands_nonhuman$pentamer_expression,0.75,na.rm = T)) & ligands_nonhuman$thymopscore > quantile(ligands_nonhuman$thymopscore,0.25,na.rm = T)] = 2
table(ligands_nonhuman$group)

#TCEM-ek, amikre nem várunk pozitivan szelektalodott T-sejteket
nptcemsd12 = unique(unname(ligands_nonhuman[ligands_nonhuman$group == 1 & ligands_nonhuman$trig == 0,"tcem"])) #
nptcemsd12_self = sapply(nptcemsd12, FUN = function(s) aLign(s, s))

#REFERENCE TCEM SET - amire varunk pozitivan szelektalodott T sejtet
load("D:/CloudStation/mygit_EXT/ImObsRecognition/objects/pentamerfreq_tcem")
load("D:/CloudStation/mygit_EXT/ImObsRecognition/objects/expr_list_tcem")
exprmeds = pbsapply(expr_list, function(x) median(x))
rm(expr_list)
exprcutoff = c(0.203075,1.095007)
load("D:/CloudStation/mygit_EXT/ImObsRecognition/objects/score_list_tpr_normalized_with_division_median_two_sides")
thymomeds = pbsapply(score_list, function(x) median(x))
rm(score_list)
thymocutoff = 0.914470
keeptcems = names(pentamerfreq_tcem)[!is.na(exprmeds) & !is.na(thymomeds)]
pentamerfreq_tcem = pentamerfreq_tcem[keeptcems]
exprmeds = exprmeds[keeptcems]
thymomeds = thymomeds[keeptcems]

imtcems = names(pentamerfreq_tcem)[pentamerfreq_tcem > 3 & exprmeds >= exprcutoff[1] & exprmeds <= exprcutoff[2] & thymomeds >= thymocutoff]
rm(exprcutoff, exprmeds, keeptcems, pentamerfreq_tcem, thymocutoff, thymomeds, ligands_nonhuman)
imtcems_self = sapply(imtcems, FUN = function(s) aLign(s, s))

no_cores = detectCores()-1
c1 = makeCluster(no_cores)
simmatrix = matrix(NA, nrow = length(imtcems), ncol = length(nptcemsd12), dimnames = list(imtcems,nptcemsd12))

for(s in nptcemsd12) {
  ind = fmatch(s,nptcemsd12)
  print(ind)
  print(Sys.time())
  clusterExport(c1, c("imtcems", "aLign", "imtcems_self", "nptcemsd12_self", "fmatch", "s"))
  simmatrix[,ind] = parSapply(cl = c1, X = imtcems, FUN = function(z) {
    aLign(s, z)/(sqrt(nptcemsd12_self[s]*imtcems_self[z]))
  })
}
save(simmatrix, file = "D:/CloudStation/mygit_EXT/ImObsRecognition/objects/simmatrix_nimtcems43")
rm(s, ind)
stopCluster(c1)
rm(c1)

View(simmatrix[1:100,1:20])
