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
This step will also denovo chimera check and filter out low abundance ZOTUs. IMPORTANT: ZOTUs with less than 3 reads will be filtered out (i.e. -minsize 3) default setting filters out ZOTUs less than 8 reads
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -unoise3 uniques_combined_merged_both_runs_fil.fa -zotus combined_merged_both_runs_zotus.fa  -tabbedout combined_merged_both_runs_zotus_report.txt -minsize 3 
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

#### 8A.1) Example job submission using [SLURM](https://wiki.hpcc.msu.edu/display/ITH/Job+Scheduling+by+SLURM)
```
nano taxa_class_OTU.sbatch

#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=20:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   taxa_class_OTU  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=EMAIL@msu.edu

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -sintax combined_merged_both_runs_otus.fa -db silva_16s_v123.fa -tabbedout combined_merged_both_runs_otus_taxonomy.sintax -strand both

###end of .sbatch

sbatch taxa_class_OTU.sbatch
```

### 8B) Classifying the ZOTUs
```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -sintax combined_merged_both_runs_zotus.fa -db silva_16s_v123.fa -tabbedout combined_merged_both_runs_zotus_taxonomy.sintax -strand both
```
This step takes approximately  19 hours and 13 minutes with ~31,000 ZOTUs against the Silva reference database

### 9) Second option for classifying taxonomy of OTUs and ZOTUs using QIIME2 against [GTDB](https://gtdb.ecogenomic.org/downloads)
It is difficult to format reference databases for use in USEARCH. QIIME2 had a Naive Bayes classifier that is much more flexible.
You lab must install QIIME2 (https://docs.qiime2.org/2019.10/install/native/) and train the classifer (https://docs.qiime2.org/2019.10/tutorials/feature-classifier/) before running the steps

#### You need to capitalize the codons before you run the classifier
```

#OTU rep set 

tr '[:lower:]' '[:upper:]'  < combined_merged_both_runs_otus.fa > combined_merged_both_runs_otus_CAP.fa
 
#ZOU rep set
 tr '[:lower:]' '[:upper:]'  < combined_merged_both_runs_zotus.fa > combined_merged_both_runs_zotus_CAP.fa

```

#### Convert these fasta files to a QIIME [artifacts](https://docs.qiime2.org/2019.10/concepts/#data-files-qiime-2-artifacts)

```
#start QIIME version which you used to train the classifier

source activate qiime2-2019.4 

#OTU
qiime tools import --type 'FeatureData[Sequence]' --input-path combined_merged_both_runs_otus_CAP.fa --output-path combined_merged_both_runs_otus.qza


#ZOTU
qiime tools import --type 'FeatureData[Sequence]' --input-path combined_merged_both_runs_zotus_CAP.fa --output-path combined_merged_both_runs_zotus.qza
```

#### Now you use the Naive Bayes classifier  
bac120_ar122_ssu_r89_ref_515R-806F_classifier.qza is the trained reference database
```
qiime feature-classifier classify-sklearn --i-classifier bac120_ar122_ssu_r89_ref_515R-806F_classifier.qza --i-reads combined_merged_both_runs_otus.qza --o-classification combined_merged_both_runs_otus_GTDBr89_taxonomy.qza 

qiime feature-classifier classify-sklearn --i-classifier bac120_ar122_ssu_r89_ref_515R-806F_classifier.qza --i-reads combined_merged_both_runs_zotus.qza --o-classification combined_merged_both_runs_zotus_GTDBr89_taxonomy.qza 
```


#### Convert the taxa table to tsv 

```
qiime tools export --input-path combined_merged_both_runs_zotus_GTDBr89_taxonomy.qza --output-path otus_GTDBr89_taxonomy

qiime tools export --input-path combined_merged_both_runs_zotus_GTDBr89_taxonomy.qza --output-path zotus_GTDBr89_taxonomy
```



### 10) Optional tree building using [PASTA](https://github.com/smirarab/pasta)
Your lab must locally install PASTA to run this step.
NOTE: This tree needs to be converted to .nwk format for use in R.
```
#modules that you need to load on HPCC
module load GNU/7.3.0-2.30  OpenMPI/3.1.1
module load Python/2.7.15

#OTU
python run_pasta.py -i combined_merged_both_runs_otus.fa -j tree_OTU -o tree_pasta_OTU/

#ZOTU
python run_pasta.py -i combined_merged_both_runs_zotus.fa -j tree_ZOTU -o tree_pasta_ZOTU/
```