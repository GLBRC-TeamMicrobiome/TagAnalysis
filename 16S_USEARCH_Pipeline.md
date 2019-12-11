# 16S USEARCH pipeline GLBRC
### Note: this pipeline merges raw fastq sequences from two illumina runs (ignore run 2 code if you’re only working with 1 run)
### Make sure you are in the folder with the extracted forward and reverse reads from Run 1

## 1) Quality checking
### 1a) First look at the quality of  raw unmerged seqs for run1
#https://www.drive5.com/usearch/manual/pipe_readprep_understand.html

`
mkdir fastq_info_run1
`

#### make a forloop to run fastx_info on every file
```
nano fasta_info_fq.sh
!#/bin/bash
for fq in *.fastq
do
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_info $fq -output fastq_info_run1/$fq
done
```

#### make file executable and run the for loop create in your `fasta_info_fq.sh` file
```
chmod +x fasta_info_fq.sh
./fasta_info_fq.sh
```

#### move to the fastq_info directory
```
cd fastq_info_run1/
```

#### Now run the below code to summarize the fastq info for all of the forward and reverse reads

```
grep "^File" * > Run1_fastq_lengths.txt
grep "^EE" * > Run1_fastq_EE.txt
```

look for any forward and reverse reads that look especially bad in terms of quality (high E is bad quality). This info will also be really helpful for troubleshooting later on (e.g. why some samples have extremely low read numbers)

### 1b) Repeat the above for Run2.



##  2) Merge the forward and reverse sequences and trim adapters (for each run individually)
### 2a) Merge Pairs
#### Make sure you are in the folder with the extracted forward and reverse reads from Run 1
https://www.drive5.com/usearch/manual/merge_options.html
#### -alnout gives you a human readable text file of the alignment and misalignments for each pair merged.
#### -tabbedout give you extensive information on the quality of the merge
#### This step takes approximately 1 minute
```
mkdir mergedfastq_run1
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_mergepairs *R1*.fastq -relabel @ -fastqout mergedfastq_run1/combined_merged_run1.fastq  -tabbedout mergedfastq_run1/combined_merged_run1_pair_report.txt -alnout mergedfastq_run1/combined_merged_run1_pair_aln.txt
```

If you are getting a low percentage of sequence pair merging (e.g. below 60%), consider trimming the reverse reads See bottom of script (“improving mergepairs”) for script on truncating reverse reads.
#### Tips on poor merging [here](https://drive5.com/usearch/manual/merge_badrev.html).

### 2b) let's check sequence quality of the merged seqs using USEARCH's [fastq_eestats2](https://www.drive5.com/usearch/manual/cmd_fastq_eestats2.html).

```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_eestats2 mergedfastq_run1/combined_merged_run1.fastq -output fastq_info/combined_merged_run1_eestats2.txt
```

### 2c) Now remove any residual bases from adapter seqs using [cut adapt](http://cutadapt.readthedocs.io/en/stable/index.html)
This code removes the forward adapter of 515F and the reverse complement of 806R. If you’re using the hpcc, load cutadapt.

```
module load cutadapt/1.8.1
cutadapt -a ATTAGAWACCCBDGTAGTCC -a GTGCCAGCMGCCGCGGTAA -o cut_combined_merged_run1.fastq combined_merged_run1.fastq > cut_adpt_results_combined_merged_run1.txt
```

#### Repeat the above merging and adapter removal for run 2.

## 3) Combine the two merged sequence files (or as many as you have Illumina runs)
### Make sure you do not use replicate sample names between the two runs
### This step is only necessary if more than one Illumina run are being analyzed together.

```
cat cut_combined_merged_run1.fastq cut_combined_merged_run2.fastq > combined_merged_both_runs.fastq
```

Before we continue, you may want to check if the sample names are formatted correctly. USEARCH does some funny cutting during the merging step. Any hyphens or underscores can be problematic and you need to remove these (use sed command and merged_cut files)

Additionally, this is a good opportunity to double check that all of your samples merged and have unique IDs using [fastx_get_sample_names](https://www.drive5.com/usearch/manual/cmd_fastx_get_sample_names.html)

```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_get_sample_names combined_merged_both_runs.fastq -output combined_merged_both_runs_samples.txt
```

## 4) Filtering and Truncate the merged seqs  to MaxEE and set length using [fastq_filter](https://www.drive5.com/usearch/manual/cmd_fastq_filter.html)

#### 250 bp is the expected overlaps with 515F and 806R
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_filter combined_merged_both_runs.fastq -fastq_maxee 1 -fastq_trunclen 250 -fastaout combined_merged_both_runs_fil.fa
```

## 5) Filter so we only have unique sequences with [fastx_uniques](https://www.drive5.com/usearch/manual/cmd_fastx_uniques.html)
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_uniques combined_merged_both_runs_fil.fa  -fastaout uniques_combined_merged_both_runs_fil.fa -sizeout
```
#### This step takes approximately 2 minutes 30 seconds with ~11 million reads

## 6) Cluster into OTUS and filter out singletons
There are two options here. **(A)** uses the traditional approach and clusters sequences into 0.97 identity cutoff OTUs. **(B)** uses unoise3 to identify ZOTUs.

### 6A) Cluster into 0.97 OTUs using UPARSE and [cluster_otus](https://www.drive5.com/usearch/manual/cmd_cluster_otus.html)
This step will also denovo chimera check and filter out singletons.You can remove single sequences prior to clustering but singletons are also removed at the OTU clustering step (cluster_otus filters out OTUs <2 and unoise3 filters ZOTUs <8)

```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -cluster_otus uniques_combined_merged_both_runs_fil.fa -otus combined_merged_both_runs_otus.fa -uparseout combined_merged_both_runs_otus_uparse.txt -relabel OTU
```

This step takes approximately 30 minutes with ~330,000 unique sequences

### 6B) Identify ZOTUs using [unoise3](https://www.drive5.com/usearch/manual/cmd_unoise3.html)
This step will also denovo chimera check and filter out low abundance ZOTUs (Less than 8 reads)
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -unoise3 uniques_combined_merged_both_runs_fil.fa -zotus combined_merged_both_runs_zotus.fa  -tabbedout combined_merged_both_runs_zotus_report.txt
```
This step takes approximately 30 minutes with ~330,000 unique sequences.
#### You must rename your representative sequence names from “Zotu” to “ZOTU” for the mapping back to the raw reads step to work correctly
```
sed -i 's/Zotu/ZOTU/g' combined_merged_both_runs_zotus.fa
```


## 7) Map reads back to OTUs at a 97% similarity score using [otutab](https://www.drive5.com/usearch/manual/cmd_otutab.html)
**-id 0.97 -strand plus are defaults**
### 7A) Mapping reads to traditional 0.97 OTUS
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -otutab combined_merged_both_runs.fastq -otus combined_merged_both_runs_otus.fa -uc combined_merged_both_runs_OTU_map.uc -otutabout combined_merged_both_runs_OTU_table.txt -biomout combined_merged_both_runs_OTU_jsn.biom -notmatchedfq combined_merged_run2_otu_unmapped.fq
```
This step takes approximately 1 hour 30 minutes with ~30,000 OTUs

### 7B) Mapping reads to ZOTUs
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -otutab combined_merged_both_runs.fastq -zotus combined_merged_both_runs_zotus.fa -uc combined_merged_both_runs_ZOTU_map.uc -otutabout combined_merged_both_runs_ZOTU_table.txt -biomout combined_merged_both_runs_ZOTU_jsn.biom -notmatchedfq combined_merged_run2_ZOTU_unmapped.fq
```
This step takes approximately 1 hour with ~31,000 ZOTUs

## 8) Classifying taxa against the reference database using [sintax](https://www.drive5.com/usearch/manual/cmd_sintax.html)
Currently used database is [silva version 123](https://www.drive5.com/usearch/manual/sintax_downloads.html).
### 8A) Classifying the traditional 0.97 OTUs
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -sintax combined_merged_both_runs_otus.fa -db silva_16s_v123.fa -tabbedout combined_merged_both_runs_otus_taxonomy.sintax -strand both
```
This step takes approximately  18 hours and 20 minutes with ~30,000 OTUs against the Silva reference database, so you may want to submit a job for it.

### 8B) Classifying the ZOTUs
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -sintax combined_merged_both_runs_zotus.fa -db silva_16s_v123.fa -tabbedout combined_merged_both_runs_zotus_taxonomy.sintax -strand both
```
This step takes approximately  19 hours and 13 minutes with ~31,000 ZOTUs against the Silva reference database
