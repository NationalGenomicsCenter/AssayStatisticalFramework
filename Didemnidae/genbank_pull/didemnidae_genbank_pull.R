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
Out_GenBank <- "didemnidae_sequences.fas" # Outputs sequence FASTA file
Out_table <- "didemnidae_records.csv" # Outputs sequence metadata (CSV)

#Gene_list <- "mitochondrion[FILT] AND complete genome[GENE]"
Gene_list <- "((coi[GENE] OR cox1 [GENE] OR cytochrome oxidase subunit I[GENE] OR cytochrome oxidase subunit 1[GENE] OR cytochrome oxidase gene[GENE]
OR coi[PRODUCT] OR cox1 [PRODUCT] OR cytochrome oxidase subunit I[PRODUCT] OR cytochrome oxidase subunit 1[PRODUCT] OR cytochrome oxidase gene[PRODUCT]
OR coi[TITL] OR cox1 [TITL] OR cytochrome oxidase subunit I[TITL] OR cytochrome oxidase subunit 1[TITL] OR cytochrome oxidase gene[TITL]) AND mitochondrion[FILT])"
#Gene_list <- "(cytb[GENE] OR cytochrome b[GENE] OR cytb[TITL] OR cytochrome b[TITL])"
#Gene_list <- "(nd1[GENE] OR nadh dehydrogenase subunit 1[GENE] OR nd1[TITL] OR nadh dehydrogenase subunit 1[TITL])"
#Gene_list <- "(nd2[GENE] OR nadh dehydrogenase subunit 2[GENE] OR nd2[TITL] OR nadh dehydrogenase subunit 2[TITL])"
#Gene_list <- "(nd4[GENE] OR nadh dehydrogenase subunit 4[GENE] OR nd4[TITL] OR nadh dehydrogenase subunit 4[TITL])"
#Gene_list <- "(atp6[GENE] OR atp synthase subunit 6[GENE] OR atp6[TITL] OR atp synthase subunit 6[TITL])"
#Gene_list <- "(d-loop[GENE] OR control region[GENE] OR control_region[GENE] OR d-loop[TITL] OR control region[TITL] OR control_region[TITL])"
#Gene_list <- "(its1[GENE] OR itsI[GENE] OR internal transcribed spacer 1[GENE] OR its1[TITL] OR itsI[TITL] OR internal transcribed spacer 1[TITL])"
#Gene_list <- "(12S[GENE] OR 12S[TITL])"
#Gene_list <- "(16S[GENE] OR 16S[TITL])"
#Gene_list <- "(matK[GENE] OR maturase K[GENE] OR matK[TITL] OR maturase K[TITL])"
#Gene_list <- "(trnL[GENE] OR trnF[GENE] OR trnL[TITL] OR trnF[TITL])"
#Gene_list <- "(mcp[GENE] OR major capsid protein[GENE] OR mcp[TITL] OR major capsid protein[TITL])"

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
  
  if(R_search$count > 0){
    
    print(paste(R_search$count, "records for", Taxa_list[i])) # Status
    
    R_fetched <- entrez_summary(db="nucleotide", 
                                web_history=R_search$web_history,
                                retmode="xml") 
    
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
  
  else{
    print(paste("No records for", Taxa_list[i])) # Status
  }
  
  print(paste(i, "of", length(Taxa_list), "done!")) # Status
}


write.csv(All_table, Out_table, row.names = FALSE) 
write(All_fasta, Out_GenBank) 

### Some summary information
All_table
species <- levels(as.factor(All_table$Species))
length(species)