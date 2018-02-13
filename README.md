# TagAnalysis
Repo for 16S and ITS workflows




-------------------------------------
# Pipeline for ITS amplicon sequences
It will be made of different scripts for each of the main bioinformatic steps, which are the following:

1) raw reads quality check and prefiltering;
- FastQC
- Bowtie2 (removing Phix reads)
- python scripts

2) read merging;
- USEARCH/VSEARCH
- FastQC

3) quality filtering and trimming;
- USEARCH
- FastQC

4) Clustering OTUs and Generating ESV (Exact Sequence Variants);
- USEARCH (UPARSE and UNOISE3 algorithms)

5) taxonomy assignments;
- CONSTAX (we developed this tool that generates a consensus taxonomy for ITS amplicon reads, please have a look https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1952-x

- RDP/SINTAX? I am happier with RDP and makes more sense to me since we have the RDP folks on campus.

# Folders
- python_scrips
- phix_index
- Reference_DB





