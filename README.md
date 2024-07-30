# Yin_CHIKV_fidelity_supp
Supplemental information for Figures 1, 4 and 5 of Mutations in Chikungunya Virus nsP4 Decrease Sensitivity  to the Broad-Spectrum Antiviral 4′-Fluorouridine.

NCBI_BioProject_info contains the metadata for BioProject PRJNA1141983.

See https://github.com/greninger-lab/RAVA_Pipeline for variant analysis used to generate figures 1 and 4.

fastp_reports contains html files showing quality filtering stats for the merged paired-end reads.

figure1 contains the following data:
 1. rava output files for the resistance profiling experiment

figure4 contains the following data:
  1. rava output files for the competition assay
     
figure5 contains the following data:

  1. Variant analysis data in vcf and csv format (in sub-directories fidelity_1 and fidelity_2), and the python script (parse_readcounts_mutrate.py) used to parse the vcfs into csv tables. 

  2. R script used to calculate mutations per 10k bases, mutated site frequency, and shannon entropy         
  (fidelity_figures_combo_final.R), and the output csv tables from that script (in sub-directory output_tables). 
