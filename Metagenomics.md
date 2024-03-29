# Biol 217 Practice Session Day -1

Login Data: ssh -X sunam236@caucluster.rz.uni-kiel.de

What we have learned so far?

1. Basic Linux 
2. Bioinformatics basic understandings 
3. Linux commands 

Block of code 

```sh
cp source destination 
```

this is the command `cp`

### Task: Practice how to upload images and links 

# Day 2 From raw reads to MAGs

## 1. Quality control of the reads

Anvio_slurm.txt - file to run the Batch Skript 

FastQC to evaluate quality of the reads
```sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc.out
#SBATCH --error=fastqc.err
#SBATCH --partition=base
#SBATCH --reservation=biol217
# Activate conda environment: 
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd /work_beegfs/sunam236/Metagenomics/0_raw_reads


for i in *fastq.gz; do fastqc $i -o  ../1_fastqc/; done
```

Output path is: Metagenomics/1_fastqc

To run the text `sbatch file.txt`
To check process `squeue -u sunam236`


Copy and open in browser to see plots:

```
scp sunam236@caucluster.rz.uni-kiel.de:/work_beegfs/sunam236/Metagenomics/1_fastqc/*.html .
```

### fastp for trimming the reads

> `--html` creates an .html report file in html format\
>`-i` R1 
>`-I` R2 
>`-R` report title, here ‘_report’ is added to each file\
>`-o` output_folder/R1.fastq.gz output file\
>`-O` output_folder/R2.fastq.gz output file\
>`-t` trim tail 1, default is 0, here 6 bases are trimmed\
>`-q` 20 reads with a phred score of <=20 are trimmed


Sample 1 
```
cd /work_beegfs/sunam236/Metagenomics/0_raw_reads

fastp -i BGR_130305_mapped_R1.fastq.gz -I BGR_130305_mapped_R2.fastq.gz -R fastp_Report -o ../2_fast/BGR_130305_mapped_R1_clean.fastq.gz -O ../2_fast/BGR_130305_mapped_R2_clean.fastq.gz -t 6 -q 20
```
Sample 2 
```
fastp -i BGR_130527_mapped_R1.fastq.gz -I BGR_130527_mapped_R2.fastq.gz -R fastp_Report -o ../2_fast/BGR_130527_mapped_R1_clean.fastq.gz -O ../2_fast/BGR_130527_mapped_R2_clean.fastq.gz -t 6 -q 20
```
Sample 3 
```
fastp -i BGR_130708_mapped_R1.fastq.gz -I BGR_130708_mapped_R2.fastq.gz -R fastp_Report -o ../2_fast/BGR_130708_mapped_R1_clean.fastq.gz -O ../2_fast/BGR_130708_mapped_R2_clean.fastq.gz -t 6 -q 20
```
fastp 
Loop: 
```
for i in `ls *_R1.fastq.gz`;
do
    second=`echo ${i} | sed 's/_R1/_R2/g'`
    fastp -i ${i} -I ${second} -R "${i}"_report -o output_folder/"${i}" -O output_folder/"${second}" -t 6 -q 20

done
```
- Fastp report shows how many reads were trimmed 

## 2. Assembly 

### megahit

- ultra-fast and memory-efficient NGS assembler optimized for metagenomes

```
cd /work_beegfs/sunam236/Metagenomics/2_fastp

megahit -1 BGR_130305_mapped_R1_clean.fastq.gz -1 BGR_130527_mapped_R1_clean.fastq.gz -1 BGR_130708_mapped_R1_clean.fastq.gz -2 BGR_130305_mapped_R2_clean.fastq.gz -2 BGR_130527_mapped_R2_clean.fastq.gz -2 BGR_130708_mapped_R2_clean.fastq.gz
--min-contig-len 1000 --presets meta-large -m 0.85 -o ../3_coassembly -t 12 

```

To view the error file whilst running: `tail file.err`


Bandage File: 

- Many small unconnected contigs
- All of them are linear, no circular contigs
- No complex structures
- The line length is proportional to the contig length

<img src="https://github.com/MarieRaasch/biol217/blob/main/Ressources/graph%20(1).png" alt="Description of the image" width="500" height="400">


### Day 3


Checking if assembly worked: 

`$ tail final.contigs.fa ` 

To count the number of contigs in your final assembly file:

`grep -c ">" final.contigs.fa`
55838

To visualize contig graph in Bandage, the first step is to convert the fasta file(s) intermediate_contigs/k{kmer_size}.contigs.fa into SPAdes-like FASTG format. The following code shows the translation from NAME.contigs.fa into NAME.fastg.

`megahit_toolkit contig2fastg 99 final.contigs.fa > final.contigs.fastg`


## Quality Assessment

**QU**ality **AS**sessment **T**ool QUAST to evaluate genome assembly.

metaquast -t 6 -o ? -m 1000 ?

Batch Script 
```
cd /work_beegfs/sunam236/Metagenomics/3_coassembly

metaquast -t 6 -o ../3_metaquast_out -m 1000
```
Looking at our report.pdf in local Dektop
```sh
scp sunam236@caucluster.rz.uni-kiel.de:/work_beegfs/sunam236/Metagenomics/3_metaquast_out/report.pdf . 
```
![Screenshot 2024-02-13 160804](https://github.com/MarieRaasch/biol217/assets/157317805/51f63da7-a11f-4f1c-aec2-7442012b4664)

![Screenshot 2024-02-13 160824](https://github.com/MarieRaasch/biol217/assets/157317805/60bb88df-6097-4dcf-bff3-dd4924613b9c)


**What is your N50 value?**     
3014

**Why is this value relevant?**

N50 is the length for which the collection of all contigs of that length or longer covers at least half an assembly.

**How many contigs are assembled?**

55838

**What is the total length of the contigs?**

142642670

## Genomes Binning 

1.  format fasta sequence IDs -> Anvi’o (which you will use later) needs this step to work properly. Run this with your contigs and with your clean reads (Hint: fastp output). As the names get changed, you need to run it on your assembly file, otherwise the names of the contigs won't match the names of the initial reads (essential for the mapping step below).

Batch Script

```sh
cd /work_beegfs/sunam236/Metagenomics/3_coassembly

anvi-script-reformat-fasta final.contigs.fa -o ../3_binning_out/contigs.anvio.fa  --min-len 1000 --simplify-names --report-file name_conversion.txt
```

## Mapping 

Then you need to map your raw reads onto your assembled contigs. Mapping will be done using bowtie2. Use the following command to index your mapping reference fasta file. Needed for the next steps and makes mapping faser.

```
cd /work_beegfs/sunam236/Metagenomics/3_binning_out

module load bowtie2
bowtie2-build contigs.anvio.fa contigs.anvio.fa.index


bowtie2 --very-fast -x contigs.anvio.fa.index -1 ../2_fastp/BGR_130305_mapped_R1_clean.fastq.gz -2 ../2_fastp/BGR_130305_mapped_R2_clean.fastq.gz -S BGR_130305.sam


bowtie2 --very-fast -x contigs.anvio.fa.index -1 ../2_fastp/BGR_130527_mapped_R1_clean.fastq.gz -2 ../2_fastp/BGR_130527_mapped_R2_clean.fastq.gz -S BGR_130527.sam


bowtie2 --very-fast -x contigs.anvio.fa.index -1 ../2_fastp/BGR_130708_mapped_R1_clean.fastq.gz -2 ../2_fastp/BGR_130708_mapped_R2_clean.fastq.gz -S BGR_130708.sam
```

--very-fast bowtie runs in very fast but less accurate end-to-end mode
-x index files with the contigs from the step before, give it the prefix name of the files (the part that comes before the dot)

-1 R1 fasta file containing the raw reads after fastp processing
-2 R2 fasta file containing the raw reads after fastp processing
-S name of the output file, don't forget the .sam part!

The output will be a sequence mapping file (SAM) with the .sam extension and which we convert to binary alignment and map (BAM) file with the .bam extension using samtools with the following loop:

SAMtools

```
module load samtools
samtools view -bS BGR_130305.sam > BGR_130305_bam_file.bam

samtools view -bS BGR_130527.sam > BGR_130527_bam_file.bam

samtools view -bS BGR_130708.sam > BGR_130708_bam_file.bam
```

# Contigs Data preparation 

"A contigs database is an anvi’o database that contains key information associated with your sequences" (anvio.org) 
"Contigs database is a modern, more talented form of a FASTA file, where you can store additional information about your sequences " 

```
cd /work_beegfs/sunam236/Metagenomics/3_binning_out
anvi-gen-contigs-database -f contigs.anvio.fa -o ../5_anvio_profiles/contigs.db -n 'biol217'
```

-f contig.fa files, used as input
-o will give you a .db file as output

When you run this command, anvi-gen-contigs-database will (documentation anvi-gen-contigs-database)

Then you need to perform an HMM search on your contigs. "Basically, in anvi’o, Hidden Markov Models (or HMMs for short) are used to search for specific genes with known functions in a larger dataset" (documentation [anvi-run-hmms] (https://anvio.org/help/7/programs/anvi-run-hmms/)

```cd /work_beegfs/sunam236/Metagenomics/5_anvio_profiles

anvi-run-hmms -c contigs.db --num-threads 4
```

Once you have your contigs database ready, and optionally your HMMs are run, you can take a quick look at it using the program anvi-display-contigs-stats:

```sh
srun --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --partition=base /bin/bash
```

```sh
conda activate anvio-8

anvi-display-contigs-stats contigs.db
```

New Terminal window: 

```sh
ssh -L 8060:localhost:8080 sunam236@caucluster.rz.uni-kiel.de
ssh -L 8080:localhost:8080 n100
```

In Browser

```http://127.0.0.1:8060```


# Binning with Anvio

## Sort and Index bam files 

```for i in *.bam; do anvi-init-bam $i -o "$i".sorted.bam; done```

## Creating anvio profile 
```

anvi-profile -i BGR_130305_bam_file.bam.sorted.bam -c ../5_anvio_profiles/contigs.db --output-dir ../5_anvio_profiles/130305

anvi-profile -i BGR_130527_bam_file.bam.sorted.bam -c ../5_anvio_profiles/contigs.db --output-dir ../5_anvio_profiles/130527

anvi-profile -i BGR_130708_bam_file.bam.sorted.bam -c ../5_anvio_profiles/contigs.db --output-dir ../5_anvio_profiles/130708

```
Merging the profiles coming from your different samples into one profile:

```
anvi-merge ./5_anvio_profiles/130305/PROFILE.db ./5_anvio_profiles/130527/PROFILE.db ./5_anvio_profiles/130708/PROFILE.db -o ./6_anvimerge -c ./5_anvio_profiles/contigs.db  --enforce-hierarchical-clustering
```
### Here you are going to use two binners Metabat2 and MaxBin2.

## Binning with Metabat2

```
anvi-cluster-contigs -p ./6_anvimerge/PROFILE.db -c ./5_anvio_profiles/contigs.db -C METABAT --driver metabat2 --just-do-it --log-file log-metabat2
anvi-summarize -p ./6_anvimerge/PROFILE.db -c ./5_anvio_profiles/contigs.db -o SUMMARY_METABAT -C METABAT
```

## Binning with Max Bin2 

```
anvi-cluster-contigs -p ./6_anvimerge/PROFILE.db -c ./5_anvio_profiles/contigs.db -C MAXBIN2 --driver maxbin2 --just-do-it --log-file log-maxbin2
anvi-summarize -p ./6_anvimerge/PROFILE.db -c ./5_anvio_profiles/contigs.db -o SUMMARY_MAXBIN2 -C MAXBIN2
```
## MAGs Quality Estimation 

Estimate your genomes completeness and contamination levels.
You can assess the quality of your bins by using

```
anvi-estimate-genome-completeness -c ./5_anvio_profiles/contigs.db -p ./6_anvimerge/PROFILE.db -C METABAT
```
If you want to check what collections you generated you can use:

```
anvi-estimate-genome-completeness -p ./6_anvimerge/PROFILE.db -c ./5_anvio_profiles/contigs.db --list-collections 
```
Visualize results:

```

anvi-interactive -p ./6_anvimerge/PROFILE.db -c ./5_anvio_profiles/contigs.db -C METABAT

```

![Metabat Bins](https://github.com/MarieRaasch/biol217/assets/157317805/ff83ee88-4bc4-4c46-9cc5-2d75150ea9e6)

## Results

Number of Archea with Metabat2: 1
1x 98,68% Coverage and 3,94% Redundancy

Number of Archea with Maxbin2: 1
94,7 % Coverage and 73,68% Redundancy 

Which binning strategy gives you the best quality for the Archaea
bins??

METABAT 98,68% Coverage and 3,94% Redundancy

Maxbin2: 94,7 % Coverage and 73,68% Redundancy 
With Metabat:
How many Archaea bins do you get that are of High Quality? 1 bin

How many Bacteria bins do you get that are of High Quality? 14 bins


# Day 4 

## Bin Refinement 

In Terminal
```
anvi-summarize -p./6_anvimerge/PROFILE.db -c ./5_anvio_profiles/contigs.db --list-collections

anvi-summarize -c ./5_anvio_profiles/contigs.db -p ./6_anvimerge/PROFILE.db  -C METABAT -o SUMMARY_METABAT2 --just-do-it

```

```
cd ./SUMMARY_METABAT2/bin_by_bin/

mkdir ../../ARCHAEA_BIN_REFINEMENT

cp ./METABAT__6-contigs.fa ../../ARCHAEA_BIN_REFINEMENT/

cp ./METABAT__10/*.fa ../../ARCHAEA_BIN_REFINEMENT/
```
```
anvi-estimate-genome-completeness -c ./5_anvio_profiles/contigs.db -p ./6_anvimerge/PROFILE.db -C METABAT > METABAT_table.txt
```

-> METABAT__6 & 10 


## METABAT 10
![Screenshot 2024-01-26 at 13-54-52 Refining METABAT__10 from METABAT ](https://github.com/MarieRaasch/biol217/assets/157317805/666d7a73-dc28-48be-954b-8810f264b413)

## METABAT 6 
![Screenshot 2024-01-26 at 13-49-20 Refining METABAT__6 from METABAT ](https://github.com/MarieRaasch/biol217/assets/157317805/bf54e95f-bf67-4770-862f-bff4f3a818b0)

# Chimera detection in MAGs

Use GUNC to check run chimera detection.

```

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate gunc

cd /work_beegfs/sunam236/Metagenomics/ARCHAEA_BIN_REFINEMENT

mkdir GUNC

for i in *.fa; do gunc run -i "$i" -r /work_beegfs/sunam236/Databases/gunc_db_progenomes2.1.dmnd --out_dir GUNC --threads 10 --detailed_output; done

gunc plot -d ./GUNC/diamond_output/METABAT__10-contigs.diamond.progenomes_2.1.out -g ./GUNC/genes_calls/gene_counts.json

gunc plot -d ./GUNC/diamond_output/METABAT__6-contigs.diamond.progenomes_2.1.out -g ./GUNC/genes_calls/gene_counts.json

```
Do you get bins that are chimeric?
hint: look at the CSS score (explained in the lecture) and the column PASS GUNC in the tables outputs per bin in your gunc_output folder.

In your own words (2 sentences max), explain what is a chimeric bin: 
- A chimeric bin, is a bin where contigs from different organisms were wrongly grouped together during binning. 

<img src="https://github.com/MarieRaasch/biol217/blob/main/Ressources/GUNC_plot.png"> 


## Manual Bin refinement 


```
cd /work_beegfs/sunam236/Metagenomics/ARCHAEA_BIN_REFINEMENT

anvi-refine -c ../5_anvio_profiles/contigs.db -C METABAT -p ../6_anvimerge/PROFILE.db --bin-id Bin_METABAT__10

anvi-refine -c ../5_anvio_profiles/contigs.db -C METABAT -p ../6_anvimerge/PROFILE.db --bin-id Bin_METABAT__6
```
Runnen im Terminal 

```
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

anvi-refine -c ../5_anvio_profiles/contigs.db -C METABAT -p ../6_anvimerge/PROFILE.db --bin-id METABAT__6

anvi-refine -c ../5_anvio_profiles/contigs.db -C METABAT -p ../6_anvimerge/PROFILE.db --bin-id METABAT__10

```
METABAT 6 looks good -> Continue with that

METABAT 10 Coverage too low -> Bin excluded for further work

how abundant are the archaea bins in the 3 samples? (relative abundance)

```anvi-inspect -p ../6_anvimerge/PROFILE.db -c ../5_anvio_profiles/contigs.db --split-nam ```

```anvi-script-get-coverage-from-bam -b  -C collection-txt -m bin```


Mean coverage: `cd /work_beegfs/sunam236/Metagenomics/ARCHAEA_BIN_REFINEMENT/bin_by_bin/METABAT__6-mean_coverage.txt`

`cat METABAT__6-mean_coverage.txt`

METABAT 6: 

BGR_130305_bam_file_bam_sorted: 8.200437007296363

BGR_130527_bam_file_bam_sorted: 5.316715747262296

BGR_130708_bam_file_bam_sorted: 3.4897028928552145

Abundance: 

![Screenshot 2024-02-14 145835](https://github.com/MarieRaasch/biol217/assets/157317805/27c942f0-1ec2-4417-bc2c-7499c07bb81b)


## Day 5 

Taxonomic assignment 

anvi-run-scg-taxonomy associates the single-copy core genes in your contigs-db with taxnomy information

```anvi-run-scg-taxonomy -c ./5_anvio_profiles/contigs.db -T 20 -P 2```

Now you can run anvi-estimate-scg-taxonomy, ‘This program makes quick taxonomy estimates for genomes, metagenomes, or bins stored in your contigs-db using single-copy core genes. 

To estimate abundance of Ribosomal RNAs within your dataset (coverage) use:

```anvi-estimate-scg-taxonomy -c ./5_anvio_profiles/contigs.db -p ../6_anvimerge/PROFILE.db --metagenome-mode --compute-scg-coverages --update-profile-db-with-taxonomy > temp.txt```

ONE final summary to get comprehensive info about your METABAT2 bins:

```anvi-summarize -p ./6_anvimerge/PROFILE.db -c ./5_anvio_profiles/contigs.db --metagenome-mode -o ./SUMMARY_METABAT2 -C METABAT2```

Did you get a species assignment to the bins previously identified?
- METABAT 6: Methanoculleus sp012797575
- METABAT 10: Methanosarcina flavescens
Does the HIGH-QUALITY assignment of the bin need revision?
According to Bowers et. al, 2017, a high quality draft has a completeness of >90 % and a contamination of <5 %, which was the case for the bin METABAT 6. So it does not need any further revision.

# Day 6 Genomics

## Assembling and annotating a Genome using Hybrid assembly

`cd $WORK/genomics`

creating output directory: 

```
mkdir 1_short_reads_qc
```
## 1. Quality control of the short reads and trimming using fastp and fastqc

```#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=5:00:00
#SBATCH --job-name=01_fastqc
#SBATCH --output=01_fastqc.out
#SBATCH --error=01_fastqc.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load micromamba/1.4.2
micromamba activate 01_short_reads_qc


# create a new folder for the output of qc 
mkdir -p $WORK/genomics/1_short_reads_qc/1_fastqc_raw
for i in $WORK/genomics/0_raw_reads/short_reads/*.gz; do fastqc $i -o $WORK/genomics/1_short_reads_qc/1_fast_pc_raw -t 32; done

jobinfo

```

```
1.2 fastp 
#mkdir -p $WORK/genomics/1_short_reads_qc/2_cleaned_reads
#fastp -i $WORK/genomics/0_raw_reads/short_reads/241155E_R1.fastq.gz \
 #-I $WORK/genomics/0_raw_reads/short_reads/241155E_R2.fastq.gz \
 #-R $WORK/genomics/1_short_reads_qc/2_cleaned_reads/fastp_report \
 #-h $WORK/genomics/1_short_reads_qc/2_cleaned_reads/report.html \
 #-o $WORK/genomics/1_short_reads_qc/2_cleaned_reads/241155E_R1_clean.fastq.gz \
 #-O $WORK/genomics/1_short_reads_qc/2_cleaned_reads/241155E_R2_clean.fastq.gz -t 6 -q 25
```
```
# create a new folder for output of qc 
mkdir -p $WORK/genomics/1_short_reads_qc/3_fastqc_cleaned
for i in $WORK/genomics/1_short_reads_qc/2_cleaned_reads/*.gz; do fastqc $i -o $WORK/genomics/1_short_reads_qc/3_fastqc_cleaned -t 32; done
```

How Good is the read quality?

Good 

How many reads did you have before trimming and how many do you have now?

R1
Before trimming: 1639549
After trimming: 1613392

R2
Before trimming: 1639549
After trimming: 1613392

Did the quality of the reads improve after trimming?

Yes

# Long reads QUality control and trimming 

## NanoPlot (QC)  & Filtlong (Trimming)

```#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=5:00:00
#SBATCH --job-name=02_long_reads_qc
#SBATCH --output=02_long_reads_qc.out
#SBATCH --error=02_long_reads_qc.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
module load micromamba/1.4.2

echo "---------long reads cleaning started---------"
eval "$(micromamba shell hook --shell=bash)"
micromamba activate 02_long_reads_qc

## 2.1 Nanoplot raw
cd $WORK/genomics/0_raw_reads/long_reads/
mkdir -p $WORK/genomics/2_long_reads_qc/1_nanoplot_raw
NanoPlot --fastq $WORK/genomics/0_raw_reads/long_reads/*.gz \
 -o $WORK/genomics/2_long_reads_qc/1_nanoplot_raw -t 32 \
 --maxlength 40000 --minlength 1000 --plots kde --format png \
 --N50 --dpi 300 --store --raw --tsv_stats --info_in_report

## 2.2 Filtlong
mkdir -p $WORK/genomics/2_long_reads_qc/2_cleaned_reads
filtlong --min_length 1000 --keep_percent 90 $WORK/genomics/0_raw_reads/long_reads/*.gz | gzip > $WORK/genomics/2_long_reads_qc/2_cleaned_reads/241155E_cleaned_filtlong.fastq.gz

## 2.3 Nanoplot cleaned
cd $WORK/genomics/2_long_reads_qc/2_cleaned_reads
mkdir -p $WORK/genomics/2_long_reads_qc/3_nanoplot_cleaned
NanoPlot --fastq $WORK/genomics/2_long_reads_qc/2_cleaned_reads/*.gz \
 -o $WORK/genomics/2_long_reads_qc/3_nanoplot_cleaned -t 32 \
 --maxlength 40000 --minlength 1000 --plots kde --format png \
 --N50 --dpi 300 --store --raw --tsv_stats --info_in_report

micromamba deactivate
echo "---------long reads cleaning completed Successfully---------"

module purge
jobinfo
```

How Good is the long reads quality?

<img src="https://https://github.com/MarieRaasch/biol217/blob/main/Ressources/Nanoplot%20before.png"> 


<img src="https://https://https://github.com/MarieRaasch/biol217/blob/main/Ressources/Nanoplot%20after.png" width="300">

How many reads do you had before trimming and how many do you have now?
Before: 15963
After: 12446

## Assembly, Quality check and Genome Annotation:
- Assembly using Quast
- Assembly quality check using Quast, CHeckM and CHeckM2
- Annotation using Prokka
- Classification using GTDBTK

```
# 3 Assembly (1 hour)-----------------------------------------------------------
echo "---------Unicycler Assembly pipeline started---------"
micromamba activate 03_unicycler
cd $WORK/genomics
mkdir -p $WORK/genomics/3_hybrid_assembly
unicycler -1 $WORK/genomics/1_short_reads_qc/2_cleaned_reads/241155E_R1_clean.fastq.gz -2 $WORK/genomics/1_short_reads_qc/2_cleaned_reads/241155E_R2_clean.fastq.gz -l $WORK/genomics/2_long_reads_qc/2_cleaned_reads/241155E_cleaned_filtlong.fastq.gz -o $WORK/genomics/3_hybrid_assembly/ -t 32
micromamba deactivate
echo "---------Unicycler Assembly pipeline Completed Successfully---------"

# 4 Assembly quality-----------------------------------------------------------
echo "---------Assembly Quality Check Started---------"

## 4.1 Quast (5 minutes)
micromamba activate 04_checkm_quast
cd $WORK/genomics/3_hybrid_assembly
mkdir -p $WORK/genomics/3_hybrid_assembly/quast
quast.py $WORK/genomics/3_hybrid_assembly/assembly.fasta --circos -L --conserved-genes-finding --rna-finding \
 --glimmer --use-all-alignments --report-all-metrics -o $WORK/genomics/3_hybrid_assembly/quast -t 32
micromamba deactivate

## 4.2 CheckM
micromamba activate 04_checkm_quast
cd $WORK/genomics/3_hybrid_assembly
mkdir -p $WORK/genomics/3_hybrid_assembly/checkm
checkm lineage_wf $WORK/genomics/3_hybrid_assembly/ $WORK/genomics/3_hybrid_assembly/checkm -x fasta --tab_table --file $WORK/genomics/3_hybrid_assembly/checkm/checkm_results -r -t 32
checkm tree_qa $WORK/genomics/3_hybrid_assembly/checkm
checkm qa $WORK/genomics/3_hybrid_assembly/checkm/lineage.ms $WORK/genomics/3_hybrid_assembly/checkm/ -o 1 > $WORK/genomics/3_hybrid_assembly/checkm/Final_table_01.csv
checkm qa $WORK/genomics/3_hybrid_assembly/checkm/lineage.ms $WORK/genomics/3_hybrid_assembly/checkm/ -o 2 > $WORK/genomics/3_hybrid_assembly/checkm/final_table_checkm.csv
micromamba deactivate

# 4.3 Checkm2
# (can not work, maybe due to insufficient memory usage)
micromamba activate 05_checkm2
cd $WORK/genomics/3_hybrid_assembly
mkdir -p $WORK/genomics/3_hybrid_assembly/checkm2
checkm2 predict --threads 32 --input $WORK/genomics/3_hybrid_assembly/* --output-directory $WORK/genomics/3_hybrid_assembly/checkm2 
micromamba deactivate
echo "---------Assembly Quality Check Completed Successfully---------"

# 5 Annotate-----------------------------------------------------------
echo "---------Prokka Genome Annotation Started---------"

micromamba activate 06_prokka
cd $WORK/genomics/3_hybrid_assembly
# Prokka creates the output dir on its own
prokka $WORK/genomics/3_hybrid_assembly/assembly.fasta --outdir $WORK/genomics/4_annotated_genome --kingdom Bacteria --addgenes --cpus 32
micromamba deactivate
echo "---------Prokka Genome Annotation Completed Successfully---------"


# 6 Classification-----------------------------------------------------------
echo "---------GTDB Classification Started---------"
# (can not work, maybe due to insufficient memory usage increase the ram in bash script)
micromamba activate 07_gtdbtk
conda env config vars set GTDBTK_DATA_PATH="$WORK/Databases/GTDBTK_day6";
micromamba activate 07_gtdbtk
cd $WORK/genomics/4_annotated_genome
mkdir -p $WORK/genomics/5_gtdb_classification
echo "---------GTDB Classification will run now---------"
gtdbtk classify_wf --cpus 12 --genome_dir $WORK/genomics/4_annotated_genome/ --out_dir $WORK/genomics/5_gtdb_classification --extension .fna 
# reduce cpu and increase the ram in bash script in order to have best performance
micromamba deactivate
echo "---------GTDB Classification Completed Successfully---------"

# 7 multiqc-----------------------------------------------------------
echo "---------Multiqc Started---------"
micromamba activate 01_short_reads_qc
multiqc -d $WORK/genomics/ -o $WORK/genomics/6_multiqc
micromamba deactivate
echo "---------Multiqc Completed Successfully---------"


module purge
jobinfo
```

Visualize the assembly using Bandage
Open Bandage and load the assembly file assembly.gfa from the assembly directory $WORK/genomics/3_hybrid_assembly/006_final_clean.gfa

![Bandage_Graph](https://github.com/MarieRaasch/biol217/assets/157317805/ef91d783-aa36-492c-82c3-161e2eb805e7)


# Multi QC report

### How good is the quality of genome?

Good read quality (98,4% passed the filter) 

### Why did we use Hybrid assembler?
To combine short reads and long reads to get a better assembly. 

What is the difference between short and long reads?

Short reads are produced by Illumina and and are about 150 bp long

### Did we use Single or Paired end reads? Why?

For long reads you do not use paired ends as it is not neccesary.

### Write down about the classification of genome we have used here

d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__Bacteroides sp002491635
				
