setwd("D:/CloudStation/mygit/ImObsRecognition/")
source("D:/CloudStation/mygit_EXT/ImObsRecognition/sequence_similarity_func.R")
create_db_proteome(blast_folder = "C:/NCBI/blast-2.11.0+/bin/", 
                   proteome_fasta_loc = "D:/CloudStation/mygit_EXT/ImObsRecognition/objects/uniprot-proteome_UP000005640+reviewed_yes.fasta", 
                   peptide_length = 9, 
                   db_name = "human_proteome_9mers")
blastp(blast_folder = "C:/NCBI/blast-2.11.0+/bin/", 
       query_peptides = c(dataset_1$epitope[nchar(dataset_1$epitope)==9],dataset_2$epitope[nchar(dataset_2$epitope)==9]), 
       db_loc = "C:/NCBI/blast-2.11.0+/bin/human_proteome_9mers.fasta", 
       out_loc = "D:/CloudStation/mygit_EXT/ImObsRecognition/objects/blastpout_nonamers", 
       threads = 6)
bl62score_nonamers = calcSeqSim(blastpout_loc = "D:/CloudStation/mygit_EXT/ImObsRecognition/objects/blastpout_nonamers")


create_db_proteome(blast_folder = "C:/NCBI/blast-2.11.0+/bin/", 
                   proteome_fasta_loc = "D:/CloudStation/mygit_EXT/ImObsRecognition/objects/uniprot-proteome_UP000005640+reviewed_yes.fasta", 
                   peptide_length = 9, 
                   db_name = "human_proteome_10mers")
blastp(blast_folder = "C:/NCBI/blast-2.11.0+/bin/", 
       query_peptides = c(dataset_1$epitope[nchar(dataset_1$epitope)==10],dataset_2$epitope[nchar(dataset_2$epitope)==10]), 
       db_loc = "C:/NCBI/blast-2.11.0+/bin/human_proteome_10mers.fasta", 
       out_loc = "D:/CloudStation/mygit_EXT/ImObsRecognition/objects/blastpout_dekamers", 
       threads = 6)
bl62score_dekamers = calcSeqSim(blastpout_loc = "D:/CloudStation/mygit_EXT/ImObsRecognition/objects/blastpout_dekamers")

bl62score = c(bl62score_nonamers, bl62score_dekamers)
save(bl62score, file = "D:/CloudStation/mygit_EXT/ImObsRecognition/objects/bl62score")