# Pipeline for processing multiplexed sequences aplified by ITS1f to ITS2.


## 1) Use [demux](https://docs.qiime2.org/2020.2/plugins/available/cutadapt/demux-paired/) in QIIME2 to demultiplex the dataset 
[To install QIIME2](https://docs.qiime2.org/2020.2/install/native/)

You first need to created a barcode file in a text editor 
Use QIIME format 
---
#SampleID	BarcodeSequence	LinkerPrimerSequence	Description
MMPRNT2002	TACGACTCTG	TTCTTGGTCATTTAGAGGAAGTAA	MMPRNT2002
---

Or use the QIIME2 format 

---
sample-id	barcode-sequence
#q2:types	categorical
MMPRNT2002	TACGACTCTG
---

Next, you need to rename your raw sequences before converting them to [QIIME2 artifacts](https://docs.qiime2.org/2020.2/concepts/)


```
mkdir QIIME2_demux
cd QIIME2_demux/
mkdir compressed_seqs
cd ..

cp RTSF_seqs_S1_L001_R1_001.fastq  QIIME2_demux/compressed_seqs/forward.fastq.gz
cp RTSF_seqs_S1_L001_R2_001.fastq  QIIME2_demux/compressed_seqs/reverse.fastq.gz
cp RTSF_seqs_S1_L001_I1_001.fastq  QIIME2_demux/compressed_seqs/barcodes.fastq.gz
```

Now convert the three fastq.gz to one qiime2 artifact using [qiime tools import](https://docs.qiime2.org/2020.2/tutorials/importing/)

```
source activate qiime2-2019.4

cd QIIME2_demux
qiime tools import --type EMPPairedEndSequences --input-path compressed_seqs --output-path multiplexed-seqs_paired_end.qza

conda deactivate
```

Now [demux](https://docs.qiime2.org/2020.2/plugins/available/cutadapt/demux-paired/) you sequences
I would recomend running a job for this step since it takes more than 5 hours for 27 million reads

```
nano Q2_demux_lib_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=10:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name Q2_demux_lib     # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@msu.edu

cd [Directory]

source activate qiime2-2019.4

qiime demux emp-paired --m-barcodes-file mapping_file.txt --m-barcodes-column BarcodeSequence --p-rev-comp-mapping-barcodes --p-no-golay-error-correction --i-seqs multiplexed-seqs_paired_end.qza --o-per-sample-sequences demultiplexed-seqs.qza --o-error-correction-details demux-details.qza
###End

sbatch Q2_demux_lib_fungi.sbatch

```

You can [visualize](https://docs.qiime2.org/2020.2/concepts/) the demux sequences

```
source activate qiime2-2019.4
qiime demux summarize  --i-data demultiplexed-seqs.qza  --o-visualization demux-seq_summary.qzv

conda deactivate
```

---
summary
Minimum:	87
Median:	143104.5
Mean:	140706.50520833334
Maximum:	320018
Total:	27015649
---

I use the online [QIME2 Viewing](https://view.qiime2.org/) since I cannot get the [qiime tools view](https://docs.qiime2.org/2020.2/tutorials/utilities/) to work on the HPCC


Now export the demux files in the form of .fastq.gz with individual files for each sample's forward and reverse reads

```
source activate qiime2-2019.4

qiime tools export  --input-path demultiplexed-seqs.qza  --output-path output/demultiplexed-seqs
#Exported demultiplexed-seqs.qza as SingleLanePerSamplePairedEndFastqDirFmt to directory output/demultiplexed-seqs

conda deactivate
```

Let's unzip the resulting demultiplex files

```
cd ~/QIIME2_demux/output/demultiplexed-seqs/
gunzip *.fastq.gz
```

## 2) Quality checking
First look at the [quality of  raw unmerged](https://www.drive5.com/usearch/manual/pipe_readprep_understand.html) seqs for the run

```
mkdir fastq_info
```

Make a forloop to run fastx_info on every file

```
nano fasta_info_fq.sh
#!/bin/bash --login
for fq in *.fastq
do
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_info $fq -output fastq_info/$fq
done
```

Make file executable and run the for loop create in your `fasta_info_fq.sh` file

```
chmod +x fasta_info_fq.sh
./fasta_info_fq.sh
```

Move to the fastq_info directory

```
cd fastq_info/
```

Now run the below code to summarize the fastq info for all of the forward and reverse reads


```
grep "^File" * > run1_fastq_lengths.txt
grep "^EE" * > run1_fastq_EE.txt
```

Look for any forward and reverse reads that look especially bad in terms of quality (high E is bad quality). This info will also be really helpful for troubleshooting later on (e.g. why some samples have extremely low read numbers)


##  3) Merge the forward and reverse sequences and trim adapters (for each run individually)
[-fastq_mergepairs](https://www.drive5.com/usearch/manual/merge_options.html)
Make sure you are in the folder with the extracted forward and reverse reads from Run 1
-alnout gives you a human readable text file of the alignment and misalignments for each pair merged.
-tabbedout give you extensive information on the quality of the merge
This step takes approximately 30 minutes

```

nano merg_seq_lib1_fungi.sbatch


#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=1:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name merge_lib     # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@msu.edu

cd [Directory]

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_mergepairs *_R1*.fastq -relabel @ -fastqout lib1_merged_fungi.fastq  -tabbedout lib1_fungi_pair_report.txt -alnout lib1_fungi_pair_aln.txt

###End
sbatch merg_seq_lib1_fungi.sbatch
#Submitted batch job 58239587

```


### Check sample names to make sure all of the samples are present


```

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_get_sample_names lib1_merged_fungi.fastq -output lib1_merged_fungi_samples.txt

```

## Repeat the above steps for each run you have.


## 4) Combine the merged sequence files from each run
Make sure you do not use replicate sample names between the runs
This step is only necessary if more than one Illumina run are being analyzed together.

```

cat [Directory with lib]/lib1_merged_fungi.fastq [Directory with lib]/lib2_merged_fungi.fastq [Directory with lib]/lib3_merged_fungi.fastq > combined_libs_fungi.fastq


```

Before we continue, you may want to check if the sample names are formatted correctly. USEARCH does some funny cutting during the merging step. Any hyphens or underscores can be problematic and you need to remove these (use sed command and merged_cut files)

Additionally, this is a good opportunity to double check that all of your samples merged and have unique IDs using [fastx_get_sample_names](https://www.drive5.com/usearch/manual/cmd_fastx_get_sample_names.html)




### You can run the USEARCH version of [phix removal](https://www.drive5.com/usearch/manual/cmd_filter_phix.html) here if you think you have phix from the sequencing



## 5) Check sequence quality of the merged seqs using USEARCH's [fastq_eestats2](https://www.drive5.com/usearch/manual/cmd_fastq_eestats2.html) then filter the sequences by quality and length

```
cd [Directory]
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_eestats2 combined_libs_fungi.fastq -output ccombined_libs_fungi_eestats2.txt
```


### Filtering and Truncate the merged seqs  to MaxEE and set length using [fastq_filter](https://www.drive5.com/usearch/manual/cmd_fastq_filter.html)

For this exmaple, 250 bp is within the expected amplification length and >91.7% of sequences meet this length cut off

```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastq_filter combined_libs_fungi.fastq -fastq_maxee 1 -fastq_minlen 250 -fastaout combined_libs_fungi_filtered.fa

```



## 5) Let's remove the primer bp and filter out any sequences missing the primer using [USEARCH primer remover](https://www.drive5.com/usearch/manual/cmd_fastx_trim_primer.html)

### Let's check to where your primers or adapters are in the fastq file
The combined dataset is very large so you need to create a random subset to test the frequency and location of primer sequences. 
Here I am using [-fastx_subsample](https://drive5.com/usearch/manual/cmd_fastx_subsample.html) to randomly choose 1 percent of the reads
```
mkdir rand_subsample

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_subsample combined_libs_fungi_filtered.fa -sample_pct 1 -fastqout rand_subsample/combined_libs_fungi_one_per.fa
```

Now make a fasta file with primers. I am using the ITS1f-ITS2 set here used by the [Earth Microbiome Project](http://www.earthmicrobiome.org/protocols-and-standards/its/)

```
nano EMP_ITS1-2_fungal_primers.fa

>ITS1f 
CTTGGTCATTTAGAGGAAGTAA
>ITS2r
GCTGCGTTCTTCATCGATGC

```
Now use [-search_oligodb](https://drive5.com/usearch/manual/cmd_search_oligodb.html) to to test the frequency and location of primer sequences.

```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -search_oligodb rand_subsample/combined_libs_fungi_one_per.fa -db EMP_ITS1-2_fungal_primers.fa -strand both -userout rand_subsample/combined_libs_fungi_one_per_output.txt -userfields query+target+diffs+tlo+thi+qlor+qhir

```

---
Info on the [userfields](https://www.drive5.com/usearch/manual/userfields.html)

query	 	Query sequence label.
target	 	Target sequenc label.
diffs	 	Number of differences (gaps + mismatches).
tlo	 	0-based start position of alignment in target sequence.
thi	 	0-based end position of alignment in target sequence.
qlo	 	0-based start position of alignment in query sequence.
qhi	 	0-based end position of alignment in query sequence.
---


## 6) Let's use the [USEARCH primer remover](https://www.drive5.com/usearch/manual/cmd_fastx_trim_primer.html) on the random sequences
For this example, useroutput from -search_oligodb puts primer sequences start at 7-15 bp. This is important for -width parameter though I usually use a larger value to be safe.


```
/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_trim_primer rand_subsample/combined_libs_fungi_one_per.fa -db EMP_ITS1-2_fungal_primers.fa -strand both  -maxdiffs 3 -width 20 -fastaout rand_subsample/combined_libs_fungi_one_per_primer_re.fa -tabbedout rand_subsample/combined_libs_fungi_one_per_primer_remove.txt

#I use grep to get a count of sequences
grep -c "^>" rand_subsample/combined_libs_fungi_one_per_primer_re.fa

```

### Now let's remove the primers from the entire fastq dataset
This takes ~30 minutes so you may want to submit a job.

```

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_trim_primer combined_libs_fungi_filtered.fa -db EMP_ITS1-2_fungal_primers.fa -strand both -maxdiffs 3 -width 20 -fastaout combined_libs_fungi_primer_filtered.fa -tabbedout combined_libs_fungi_filtered_primer_re_output.txt

#I use grep to get a count of sequences
grep -c "^>" combined_libs_fungi_primer_filtered.fa

```


## 7) Filter so we only have unique sequences with [fastx_uniques](https://www.drive5.com/usearch/manual/cmd_fastx_uniques.html)

```

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -fastx_uniques combined_libs_fungi_primer_filtered.fa   -fastaout uniques_combined_libs_fungi_primer_filtered.fa -sizeout

```


## 8) Cluster into OTUS and filter out singletons
There are two options here. **(A)** uses the traditional approach and clusters sequences into 0.97 identity cutoff OTUs. **(B)** uses unoise3 to identify ZOTUs.

### 8A) Cluster into 0.97 OTUs using UPARSE and [cluster_otus](https://www.drive5.com/usearch/manual/cmd_cluster_otus.html)
This step will also denovo chimera check and filter out singletons.You can remove single sequences prior to clustering but singletons are also removed at the OTU clustering step (cluster_otus filters out OTUs <2 and unoise3 filters ZOTUs <8)

This step takes approximately 30 minutes with ~330,000 unique sequences

### 8B) Identify ZOTUs using [unoise3](https://www.drive5.com/usearch/manual/cmd_unoise3.html)
This step will also denovo chimera check and filter out low abundance ZOTUs (Less than 8 reads)
This step takes approximately 30 minutes with ~330,000 unique sequences.
I am going to combine step eight and submit a job to slurm.

```

nano OTU_ZOTU_clustering_fungi.sbatch

#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=25:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   clust_OTU_ZOTU_Fung  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@msu.edu

cd [Directory]

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -cluster_otus uniques_combined_libs_fungi_primer_filtered.fa -otus rep_set_combined_libs_fungi_otus.fa -relabel OTU

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -unoise3 uniques_combined_libs_fungi_primer_filtered.fa -zotus  rep_set_combined_libs_fungi_ZOTUs.fa -minsize 3 -relabel ZOTU

###end of .sbatch

sbatch OTU_ZOTU_clustering_fungi.sbatch

```



## 9) Map reads back to OTUs at a 97% similarity score using [otutab](https://www.drive5.com/usearch/manual/cmd_otutab.html)
**-id 0.97 -strand plus are defaults**
### 9A) Mapping reads to traditional 0.97 OTUS

This step takes approximately 1 hour 30 minutes with ~30,000 OTUs

### 9B) Mapping reads to ZOTUs

This step takes approximately 1 hour with ~31,000 ZOTUs


```

nano OTU_ZOTU_Fungi_mapping.sbatch

#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########

#SBATCH --time=45:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   mapping_OTU_ZOTU  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@msu.edu

cd [Directory]

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -otutab combined_libs_fungi.fastq -otus rep_set_combined_libs_fungi_otus.fa -otutabout combined_libs_fungi_OTU_table.txt -notmatchedfq combined_libs_fungi_otu_unmapped.fq

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -otutab combined_libs_fungi.fastq -zotus rep_set_combined_libs_fungi_ZOTUs.fa -otutaboutcombined_libs_fungi_ZOTUs_table.txt -notmatchedfq combined_libs_fungi_ZOTUs_unmapped.fq



###end of .sbatch

sbatch OTU_ZOTU_Fungi_mapping.sbatch

```

## 10) Classifying taxa against the reference database using [sintax](https://www.drive5.com/usearch/manual/cmd_sintax.html)


### 10A) Classifying the traditional 0.97 OTUs
### 10B) Classifying the ZOTUs

```
nano taxa_class_MMPRNT018_Z.OTU.sbatch

#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########
 
#SBATCH --time=45:30:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=5           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=20G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   taxa_class_Z.OTU  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@msu.edu
cd [Directory]

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -sintax rep_set_combined_libs_fungi_otus.fa -db utax_reference_dataset_all_02.02.2019.fasta -tabbedout taxonomy_combined_libs_fungi_otus.sintax -strand both

/mnt/research/rdp/public/thirdParty/usearch11.0.667_i86linux64 -sintax rep_set_combined_libs_fungi_ZOTUs.fa -db utax_reference_dataset_all_02.02.2019.fasta -tabbedout taxonomy_combined_libs_fungi_ZOTUs.sintax -strand both

###end of .sbatch

sbatch taxa_class_MMPRNT018_Z.OTU.sbatch
#Submitted batch job 58960105
```

## Alternatively you can classify using CONSTAX which uses sintax, utax, and RDP to create a consensus taxonomy

[CONSTAX](https://github.com/Gian77/CONSTAX) in available in the HPCC at MSU.
You cab find CONSTAX v.2 here /mnt/research/common-data/Bio/UserDownloads/CONSTAX

You needed to create a conda environment. (Evans Lab can use source activate CONSTAX_py2 location /mnt/research/EvansLab/Software/anaconda2/envs/CONSTAX_py2/bin)

```

nano contax_class_fungi_OTU.sbatch

#!/bin/bash -login

#SBATCH --time=01:30:00	            # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=1-2                # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=2                  # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=20           # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=32G            # memory required per allocated CPU (or core) - amount of memory (in bytes)
#SBATCH --job-name   constax_OTU_Fung  # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@msu.edu


# activate your conda environment with python 2.7 and Bio-python here

source activate CONSTAX_py2

# set names for  database, otu file, and confidence level 
sed -i 's/ref_database=.*/ref_database="sh_general_release_fungi_35077_RepS_04.02.2020.fasta"/' /mnt/research/common-data/Bio/UserDownloads/CONSTAX/config

###### ERROR 1 - here just use the name of your actual otu_table, not the full PATH.
sed -i 's/otu_file=.*/otu_file="rep_set_combined_libs_fungi_otus.fa"/' /mnt/research/common-data/Bio/UserDownloads/CONSTAX/config
sed -i 's/conf=.*/conf="0.8"/' /mnt/research/common-data/Bio/UserDownloads/CONSTAX/config

# set "perform_training=yes" if you want to retrain the reference or you are using a different reference.
sed -i 's/perform_training=.*/perform_training=no/' /mnt/research/common-data/Bio/UserDownloads/CONSTAX/config

# In the latter, you must copy the new reference in the CONSTAX/DB repo. Similar code line of that below
# used for copying the otu_table into the CONSTAX/otus/ repo.

##### ERROR 2 - Here, you actually have to copy your otu_table in the otus/ folder inside the CONSTAX hpcc folder
# cp otus/otus_ITS2_R1.fasta /mnt/research/common-data/Bio/UserDownloads/CONSTAX/otus/
cp /[OTU Directory]/rep_set_combined_libs_fungi_otus.fa /mnt/research/common-data/Bio/UserDownloads/CONSTAX/otus/

# Running constax.sh here
cd /mnt/research/common-data/Bio/UserDownloads/CONSTAX/
sh /mnt/research/common-data/Bio/UserDownloads/CONSTAX/constax.sh
mv /mnt/research/common-data/Bio/UserDownloads/CONSTAX/outputs/ /[Your Directory]
rm /mnt/research/common-data/Bio/UserDownloads/CONSTAX/taxonomy_assignments/ -rf

source deactivate

sbatch contax_class_fungi_OTU.sbatch

```
