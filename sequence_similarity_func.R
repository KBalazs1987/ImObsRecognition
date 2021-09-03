#Sequence similarity function

#Create Database from a proteome
# blast_folder = "C:/NCBI/blast-2.11.0+/bin/" #directory name does not contain a white space (e.g. "Program Files" is not allowed)
# proteome_fasta_loc = "D:/CloudStation/uniprot/uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_+AND+review--.fasta"
# peptide_length = 15
# db_name = "human_proteome_15mers"
create_db_proteome <- function(blast_folder, proteome_fasta_loc, peptide_length, db_name) {
  proteome = protr::readFASTA(proteome_fasta_loc)
  mers = sapply(proteome, FUN = function(z) substring(z, 1:(nchar(z)-(peptide_length - 1)), peptide_length:nchar(z)))
  mers = unique(unlist(mers))
  strings = rbind(paste0(">", mers), mers)
  strings = as.vector(strings)
  writeLines(strings, paste0(blast_folder, db_name, ".fasta"))
  message("Makeblastdb application is curently produces BLAST databases from your proteome...")
  system2(command = paste0(blast_folder, "makeblastdb.exe"), args = c(paste0("-in ", blast_folder, db_name, ".fasta"), "-dbtype prot", "-title epitopes"))
  message(paste0("Database is created in ", blast_folder))
}

#Create databse from a peptide list
# blast_folder = "C:/NCBI/blast-2.11.0+/bin/"
# peptides: character vector of peptides
# db_name = "own_peptides"
create_db_peptides <- function(blast_folder, peptides, db_name) {
  strings = rbind(paste0(">", peptides), peptides)
  strings = as.vector(strings)
  writeLines(strings, paste0(blast_folder, db_name, ".fasta"))
  message("Makeblastdb application is curently produces BLAST databases from your peptides...")
  system2(command = paste0(blast_folder, "makeblastdb.exe"), args = c(paste0("-in ", blast_folder, db_name, ".fasta"), "-dbtype prot", "-title epitopes"))
  message(paste0("Database is created in ", blast_folder))
}

#BLASTP
# blast_folder = "C:/NCBI/blast-2.11.0+/bin/"
# query_peptides: character vector of peptides
# db_loc = "C:/NCBI/blast-2.11.0+/bin/human_proteome_15mers.fasta"
# out_loc = "C:/NCBI/blast-2.11.0+/bin/out_sample3"
# threads = 7
blastp <- function(blast_folder, query_peptides, db_loc, out_loc, threads = 1, keep_query_fasta_peptides = "NO") {
  strings = rbind(paste0(">", query_peptides), query_peptides)
  strings = as.vector(strings)
  writeLines(strings, paste0(blast_folder, "query_peptides.fasta"))
  message("blastp is already running...")
  system2(command = paste0(blast_folder, "blastp.exe"), 
          args = c(paste0("-query ", paste0(blast_folder, "query_peptides.fasta")), 
                   paste0("-db ", db_loc), 
                   paste0("-out ", out_loc),
                   paste0("-num_threads ", threads),
                   "-outfmt 6",
                   "-task blastp-short",
                   "-ungapped",
                   "-comp_based_stats F",
                   "-evalue 10000000",
                   "-max_target_seqs 100"))
  message(paste0("Sequences were compared. The location of the output file: ", out_loc))
  if(keep_query_fasta_peptides == "NO") unlink(paste0(blast_folder, "query_peptides.fasta"))
}

#Find the maximum of sequence similarity
#blastpout_loc = "C:/NCBI/blast-2.11.0+/bin/out_sample3"
calcSeqSim = function(blastpout_loc) {
  aLign = function(seq_1, seq_2, mtx = protr::AABLOSUM62, aas = colnames(protr::AABLOSUM62)) {
    sum(mtx[cbind(fmatch(unlist(strsplit(seq_1, "")), aas), fmatch(unlist(strsplit(seq_2, "")), aas))])
  }
  library(pbapply)
  library(fastmatch)
  blastpout = read.table(blastpout_loc, header = F, stringsAsFactors = F)
  mers = unique(blastpout$V1)
  ss_scores = pbsapply(mers, function(x) {
    max(sapply(blastpout$V2[blastpout$V1 == x], function(y) {
      aLign(x, y)/(sqrt((aLign(x, x))*(aLign(y, y))))
    }), na.rm = T)
  })
  return(ss_scores)
}
