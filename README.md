# Fusion_junction
Analysis of junctino sequences of fusions in colororactal cancer samples

# Terminology and notes
CRC: Colorectal Cancer
nonCRC: cancer of types other than Colorectal Cancer
MSIstatus: Microsatellite instability of tumors; can be stable (MSS), highly instable (MSI_H) or ambiguous
  * ambitous samples were exculeded from this analysis
Junction sequence: 

# Raw data
Junction_Sequences_CRC_4_19_20_deid.xlsx  # 245 Gene fusion events in CRC tumors (MSI-H, MSS, two MSI ambiguous results to be ruled out)
Junction_Sequences_NTRK1_nonCRC_4_19_20_deid.xlsx  # 112 Gene fusions events involving NTRK1 in nonCRC tumors (three MSI_H cases, 109 MSS cases)
data_Junction_Regions_CRC_4_19_20_deid.fa  # FASTA record of 101 bp sequences of genes involved in Junction 1 in CRC samples
data_Junction_Regions_CRC_4_19_20_deid_J2.fa  # FASTA record of 101 bp sequences of genes involved in Junction 2 in CRC samples
ucsc_genes.tsv  # Gene info downloaded form UCSC genome browser, hg19.


# Check GC content and dinucleotide patterns in junction sequence in CRC tumors
Fusion_Gene_CRC.ipynb

# Check GC content and dinucleotide patterns in junction sequence involving gene NTRK1 in both CRC and nonCRC tumors
Fusion_Gene_ NTRK1_partners.ipynb
Fusion_Gene_nonCRC.ipynb

# Create BED files for visualization of breakpoints in UCSC Genome browser
Bed_files_for_USCS_plot.ipynb

# Plot GC content up- and downstream of breakpoint in junction sequence
Strand_specificity_plot.ipynb
Strand_specificity_plot2.ipynb

# Plot GC content up- and downstream of breakpoint in reference gene sequence with length of 101 bp
Window_boxplot_J1.ipynb
Window_boxplot_J2.ipynb
Window_lineplot.ipynb
