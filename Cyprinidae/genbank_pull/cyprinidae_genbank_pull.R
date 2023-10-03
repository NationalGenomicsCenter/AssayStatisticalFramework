##################################################################################################
### Script outputs one CSV and one FASTA for all taxa; TXT file cannot contain > ~75 taxa
##################################################################################################
### Load packages and set up data frames
##################################################################################################

library(rentrez) 
library(stringr)
library(dplyr)

### Set working directory and specify input and output files
Taxa_list <- file.choose() # List of species names (one per row)...no parentheses or quotes allowed
Out_GenBank <- "cyprinidae_sequences.fas" # Outputs sequence FASTA file
Out_table <- "cyprinidae_records.csv" # Outputs sequence metadata (CSV)

Gene_list <- "((cytb[GENE] OR cytochrome b[GENE]
OR cytb[PRODUCT] OR cytochrome b[PRODUCT]
OR cytb[TITL] OR cox1 [TITL] OR cytochrome b[TITL] AND mitochondrion[FILT])"


##################################################################################################
### Make GenBank search query and retrieve entries
##################################################################################################
### Remove unwanted spaces from TXT file
Taxa_list <- read.table(Taxa_list, sep="\n")$V1
# Taxa_list <- gsub("?", " ", Taxa_list)
# Taxa_list <- gsub("?", "", Taxa_list)
Taxa_list <- trimws(Taxa_list)

### Empty files to populate in the loop below
All_table <- data.frame(matrix(ncol = 5, nrow = 0)) # records
colnames(All_table) <- c("Accession_title","UID","Species","Length","Location")

All_fasta <- c() # fastas

### REVISION: Loop through each of the taxa
for(i in 1:length(Taxa_list)){
  
  ### Specify search term
  Orgn_term <- paste('"', Taxa_list[i], '"[ORGN] OR ', '"', Taxa_list[i], '"[TITL] OR ', sep="") # Swap ORGN for ACCN to search for accession numbers
  
  Genbank_term <- paste(paste(Orgn_term, collapse=""),
                        "NULL", "[ORGN] AND ", Gene_list, sep="") # Use entrez_db_searchable("nucleotide") for a full list of search terms
  
  # Swap ORGN for ACCN to search for accession numbers
  # Add "100:1000[SLEN] AND " or equivalent to specify sequence length range to return
  
  ### Retrieve GenBank sequence entries
  R_search <- entrez_search(db="nucleotide",                    
                            term=Genbank_term,
                            use_history = T) 
  
  R_fetched <- entrez_summary(db="nucleotide", 
                              web_history=R_search$web_history,
                              retmode="xml") 
  
  print(paste(i, "of", length(Taxa_list), "fetched!")) # Status
  
  if(length(R_fetched) > 0){
    
    print(paste(length(R_fetched), "records for", Taxa_list[i])) # Status
    
    ### Make list of lists with records
    Accession_GenBank <- c()
    Title <- c()
    Basepairs <- c()
    UID_GenBank <- names(R_fetched)
    Orgn_GenBank <- c()
    Location_GenBank <- c()
    
    for(j in 1:length(R_fetched)){
      Accession_GenBank[j] <- R_fetched[[j]]$Caption  
      Title[j] <- R_fetched[[j]]$Title                                          
      Basepairs[j] <- R_fetched[[j]]$Slen  
      Orgn_GenBank[j] <- R_fetched[[j]]$Organism 
      ifelse(is.null(R_fetched[[j]]$SubName)==T,
             Location_GenBank[j] <- NA,
             Location_GenBank[j] <- R_fetched[[j]]$SubName)                  
    }
  }
  
  else{
    print(paste("No records for", Taxa_list[i])) # Status
  }
  
  Accession_title <- paste(Accession_GenBank, Title, sep="")
  
  ### Records file
  Records_GenBank <- data.frame(cbind(Accession_title=Accession_title,
                                      Accession=Accession_GenBank,  
                                      Title=Title,                                      
                                      Length=Basepairs,   
                                      UID=UID_GenBank,                          
                                      Organism=Orgn_GenBank,                    
                                      Location=Location_GenBank))     
  
  Records_GenBank$Species <- word(Records_GenBank$Title, 1, 2)
  Species <- Records_GenBank$Species
  Records_GenBank <- Records_GenBank[order(Species),] 
  Records_GenBank <- Records_GenBank[, c(1, 5, 8, 4, 7)]
  
  All_table <- rbind(All_table, Records_GenBank) # stick them together
  
  ### FASTA download
  IDs <- Records_GenBank$UID 
  chunk_size <- 100
  IDs_chunked <- split(IDs, ceiling(seq_along(IDs)/chunk_size))
  
  Genbank_fasta <- c()
  
  for(k in 1:length(IDs_chunked)){
    Genbank_fasta <- append(Genbank_fasta,
                            entrez_fetch(db="nucleotide",
                                         id=IDs_chunked[[k]],
                                         rettype="fasta"))
    print(paste("Chunk",k,"of",length(IDs_chunked),"for",Taxa_list[i],"done!"))
  }
  
  All_fasta <- append(All_fasta, Genbank_fasta)
}


write.csv(All_table, Out_table, row.names = FALSE) 
write(All_fasta, Out_GenBank) 
