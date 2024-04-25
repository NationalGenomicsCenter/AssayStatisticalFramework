# AssayStatistcialSpecificityInference

Analysis for Wilcox et al. THE UNKNOWN UNKNOWN: A FRAMEWORK FOR ASSESSING ENVIRONMENTAL DNA ASSAY SPECIFICITY AGAINST UNSAMPLED TAXA https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13932

We assembled species lists, sequence data, and eDNAssay assignment probability for four assays

### Nemouridae (L. tumana assay)
/genbank_pull = files for pull of sequence data from GenBank August 2023
- nemouridae_genbank_pull = R script; uses nemouridae_species_list.txt as input
- nemouridae_records = sequence information from GenBank
- nemouridae_sequences = fastas from GenBank
- nemouridae_sequences_aligned = fastas after CLUSTAL alignment

/ednassay = files for running eDNAssay model
- LEDN_eDNAssay = R script
- nemouridae_with_oligos = input alignment for eDNAssay
- LEDN_MMs = alignment information parsed by the script for the random forest model
- LEDN_APs = eDNAssay output
	
nemouridae_species_list = COL list of extant Nemouridae taxa
LEDN_summary = summarized eDNAssay assignment probabilities (max) by taxon, with sample size per taxon and genus membership

### Didemnidae (D. perlucidum assay)
/genbank_pull = files for pull of sequence data from GenBank August 2023
- didemnidae_genbank_pull = R script; uses didemnidae_species_list.txt as input
- didemnidae_records = sequence information from GenBank
- didemnidae_sequences = fastas from GenBank
- didemnidae_sequences_aligned = fastas after CLUSTAL alignment

/ednassay = files for running eDNAssay model
- DIME_eDNAssay = R script
- didemnidae_with_oligos = input alignment for eDNAssay
- DIME_APs = eDNAssay output
	
didemnidae_species_list = COL list of extant Didemnidae taxa
DIME_summary = summarized eDNAssay assignment probabilities (max) by taxon, with sample size per taxon and genus membership

### Cambaridae (F. virilis assay)
/genbank_pull = files for pull of sequence data from GenBank August 2023
- cambaridae_genbank_pull = R script; uses  cambaridae_species_list.txt as input
- cambaridae_records = sequence information from GenBank
- cambaridae_sequences = fastas from GenBank
- cambaridae_one_representative = thinned to one representative per taxon
- cambaridae_one_representative_aligned = after CLUSTAL alignment

/ednassay = files for running eDNAssay model
- cambaridae_one_representative_oligos = input alignment for eDNAssay

cambariade_species_list = COL list of extant Cambaridae taxa; note that some taxa originally listed as "Orconectes" which are now recognized as "Faxonius" were changed to reflect current taxonomy
VIR2_APs = eDNAssay assignment probabilities (one representative each taxon) 

### Cyprinidae (C. carpio assay)
/genbank_pull = files for pull of sequence data from GenBank August/Sept 2023
- cyprinidae_genbank_pull = R script; uses  cambaridae_species_list.txt as input
- cyprinidae_records = sequence information from GenBank
- cyprinidae_sequences = fastas from GenBank
- cyprinidae_one_representative = thinned to one representative per taxon

/ednassay = files for running eDNAssay model
- due to the large number of sequences, we aligned and predicted assignment probability in batches (1-5 + supplement); n = 6 files with assignment probabilities and n = 12 with alignments and oligo alignments (fasta)

cyprinidae_species_list = COL list of extant Cyprinidae taxa
CMCP_all_APs = eDNAssay assignment probabilities (one rep each taxon), pulling from files in /ednassay

### Analysis
- summaries LEDN_, DIME_, VIR2_, CMCP_summary files here; some manual edits for consistent formatting to pull into R scripts
- lednia, didemnum, cambaridae, and cyprinidae R scripts for analyses for each assay (comments in-line for scripts); pulls in script functions.r
- output files from R from four scripts above LEDN_, DIME_, VIR2_, CMCP_
- figures_script for generating figures; images pasted in manually
- appendix_2.R script for comparing weighting schemes; uses input file CMCP_genera.csv

### BiolerPlate
Files here are intended to be a user-friendly way to implement the framework described in the main manuscript. This Markdown file generates a standardized format report. Note that the report does not conduct a goodness-of-fit test to assess suitability of the beta distribution fit to the provided data or conduct rarifaction or bootstrapping of data. Both of these may be useful.
- Specificity_biolerplate = R script and example output file (PDF)
- Example_Genus_Species_AP = example input file for the Markdown boilerplate
