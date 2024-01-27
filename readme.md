# Biol 217 Practice Session 

# Day -1

Login Data: ssh -X sunam236@caucluster.rz.uni-kiel.de
Password: biol217_2024!!

Linux introduction: Basic commands in Linux Terminal: 

- pwd
- ls
- cd
- mkdir 
- cp
- mv
- cat 

Github notes: 

How to write block of code 

```sh
cp source destination 
```

this is the command `cp`


# Day-2 

### General Information: 

To run commands in Terminal: use bash scripts

Creating a file for bash Scripts: 

Anvio_slurm.txt 

Command to run bash script: `sbatch <jobscript> `

Command to manage submitted jobs: `squeue -u <username>`

Command terminate running job: `scancel <jobid>`

Viewing Error file whilst bash script job is still running: `tail file.err`

### Step 1. 
Pre-processing the raw reads (trimming adapters, removing chimeras, removing phiX174 virus sequences…)

## Quality control
Evaluate quality of sequenced data
### FastQC
- Gives Phred Quality score

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

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd /work_beegfs/sunam236/Metagenomics/0_raw_reads


for i in *fastq.gz; do fastqc $i -o  ../1_fastqc/; done
```

Output path is: Metagenomics/1_fastqc

Copy and open in web browser: 

```
scp sunam236@caucluster.rz.uni-kiel.de:/work_beegfs/sunam236/Metagenomics/1_fastqc/*.html .
```

### Fastp
- Process reads
- R1 and R2 as we have paired ends reads

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


## Assembly 
- Perform genome assembly 
### Megahit
- Uses Fastp output 
```
cd /work_beegfs/sunam236/Metagenomics/2_fastp

megahit -1 BGR_130305_mapped_R1_clean.fastq.gz -1 BGR_130527_mapped_R1_clean.fastq.gz -1 BGR_130708_mapped_R1_clean.fastq.gz -2 BGR_130305_mapped_R2_clean.fastq.gz -2 BGR_130527_mapped_R2_clean.fastq.gz -2 BGR_130708_mapped_R2_clean.fastq.gz
--min-contig-len 1000 --presets meta-large -m 0.85 -o ../3_coassembly -t 12 

```

# (Day 3)


Checking if assembly worked: 

`$ tail final.contigs.fa ` 

To count the number of contigs in your final assembly file:

`grep -c ">" final.contigs.fa`
55838

To visualize contig graph in Bandage, the first step is to convert the fasta file(s) intermediate_contigs/k{kmer_size}.contigs.fa into SPAdes-like FASTG format. The following code shows the translation from NAME.contigs.fa into NAME.fastg.

`megahit_toolkit contig2fastg 99 final.contigs.fa > final.contigs.fastg`


## Quality Assessment 

### Metaquast

- **QU**ality **AS**sessment **T**ool to evaluate genome assembly.


Bash Script 
```
cd /work_beegfs/sunam236/Metagenomics/3_coassembly

metaquast -t 6 -o ../3_metaquast_out -m 1000
```
Looking at our report.pdf in local Desktop
```sh
scp sunam236@caucluster.rz.uni-kiel.de:/work_beegfs/sunam236/Metagenomics/3_metaquast_out/report.pdf . 
```

**What is your N50 value?**     
3014

**Why is this value relevant?**

N50 is the length for which the collection of all contigs of that length or longer covers at least half an assembly.

**How many contigs are assembled?**

55838

**What is the total length of the contigs?**

142642670

# Genomes Binning 

### Reformatting data for Anvio
The first thing you will do is format your fasta sequence IDs. Anvi’o (which you will use later) needs this step to work
properly. Run this with your contigs and with your clean reads (Hint: fastp output). As the names get changed, you need to run it on your assembly file, otherwise the names of the contigs won't match the names of the initial reads (essential for the mapping step below).

Bash Script

```sh
cd /work_beegfs/sunam236/Metagenomics/3_coassembly

anvi-script-reformat-fasta final.contigs.fa -o ../3_binning_out/contigs.anvio.fa  --min-len 1000 --simplify-names --report-file name_conversion.txt
```

## Mapping 
- Mapping raw reads onto assembled contigs

### Bowtie2

```cd /work_beegfs/sunam236/Metagenomics/3_binning_out

module load bowtie2
bowtie2-build contigs.anvio.fa contigs.anvio.fa.index


bowtie2 --very-fast -x contigs.anvio.fa.index -1 ../2_fastp/BGR_130305_mapped_R1_clean.fastq.gz -2 ../2_fastp/BGR_130305_mapped_R2_clean.fastq.gz -S BGR_130305.sam


bowtie2 --very-fast -x contigs.anvio.fa.index -1 ../2_fastp/BGR_130527_mapped_R1_clean.fastq.gz -2 ../2_fastp/BGR_130527_mapped_R2_clean.fastq.gz -S BGR_130527.sam


bowtie2 --very-fast -x contigs.anvio.fa.index -1 ../2_fastp/BGR_130708_mapped_R1_clean.fastq.gz -2 ../2_fastp/BGR_130708_mapped_R2_clean.fastq.gz -S BGR_130708.sam
```


The output will be a sequence mapping file (SAM) with the .sam extension and which we convert to binary alignment and map (BAM) file with the .bam extension using samtools with the following loop

### SAMtools

```
module load samtools
samtools view -bS BGR_130305.sam > BGR_130305_bam_file.bam

samtools view -bS BGR_130527.sam > BGR_130527_bam_file.bam

samtools view -bS BGR_130708.sam > BGR_130708_bam_file.bam
```

# Contigs Data preparation 

- Compute k-mer frequencies
- Soft split contrigs longer than 20,000 bp
- Identify ORFs

```
cd /work_beegfs/sunam236/Metagenomics/3_binning_out
anvi-gen-contigs-database -f contigs.anvio.fa -o ../5_anvio_profiles/contigs.db -n 'biol217'
```

### HMM search
- Hidden Markov Models: search for specific genes with known functions in a larger dataset

```cd /work_beegfs/sunam236/Metagenomics/5_anvio_profiles

anvi-run-hmms -c contigs.db --num-threads 4
```

Display: 

```sh
srun --reservation=biol217 --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --nodelist=n100 --partition=base /bin/bash
```

```sh
conda activate anvio-8

anvi-display-contigs-stats contigs.db
```

In a new Terminal window: 

```sh
ssh -L 8060:localhost:8080 sunam236@caucluster.rz.uni-kiel.de
ssh -L 8080:localhost:8080 n100
```

In a browser tab:

```http://127.0.0.1:8060```


# Binning with Anvio

### Sort and Index bam files 

```for i in *.bam; do anvi-init-bam $i -o "$i".sorted.bam; done```

## Anvi-profile
- Recovery of mean coverage
- Characterisation of SNVs 

```

anvi-profile -i BGR_130305_bam_file.bam.sorted.bam -c ../5_anvio_profiles/contigs.db --output-dir ../5_anvio_profiles/130305

anvi-profile -i BGR_130527_bam_file.bam.sorted.bam -c ../5_anvio_profiles/contigs.db --output-dir ../5_anvio_profiles/130527

anvi-profile -i BGR_130708_bam_file.bam.sorted.bam -c ../5_anvio_profiles/contigs.db --output-dir ../5_anvio_profiles/130708

```
## Anvi-merge
Merging the profiles coming from your different samples into one profile:


```
anvi-merge ./5_anvio_profiles/130305/PROFILE.db ./5_anvio_profiles/130527/PROFILE.db ./5_anvio_profiles/130708/PROFILE.db -o ./6_anvimerge -c ./5_anvio_profiles/contigs.db  --enforce-hierarchical-clustering
```
# Binning
- Here you are going to use two binners Metabat2 and MaxBin2.

### Binning with Metabat2

```
anvi-cluster-contigs -p ./6_anvimerge/PROFILE.db -c ./5_anvio_profiles/contigs.db -C METABAT --driver metabat2 --just-do-it --log-file log-metabat2
anvi-summarize -p ./6_anvimerge/PROFILE.db -c ./5_anvio_profiles/contigs.db -o SUMMARY_METABAT -C METABAT
```

### Binning with Max Bin2 

```
anvi-cluster-contigs -p ./6_anvimerge/PROFILE.db -c ./5_anvio_profiles/contigs.db -C MAXBIN2 --driver maxbin2 --just-do-it --log-file log-maxbin2
anvi-summarize -p ./6_anvimerge/PROFILE.db -c ./5_anvio_profiles/contigs.db -o SUMMARY_MAXBIN2 -C MAXBIN2
```
## MAGs Quality Estimation 

- Estimate your genomes completeness and contamination levels.
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

## Results

Number of Archea with Metabat2: 1
1x 98,68% Coverage and 3,94% Redundancy

Number of Archea with Maxbin2: 1
94,7 % Coverage and 73,68% Redundancy 

Which binning strategy gives you the best quality for the Archaea
bins??

METABAT 98,68% Coverage and 3,94% Redundancy

Maxbin2: 94,7 % Coverage and 73,68% Redundancy 

How many Archaea bins do you get that are of High Quality? 1

1

How many Bacteria bins do you get that are of High Quality?

# Day 4 

## Bin Refinement 

- Gettig a list of our collection

In Terminal: 
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
anvi-estimate-genome-completeness -c ./5_anvio_profiles/contigs.db -p ./6_anvimerge/PROFILE.db -C METABAT > METABAT_table.txt?
```

Chosen Archaea bins: 

-> METABAT__6 & 10 

# Chimera detection in MAGs

### GUNC

- Detection of chimerisms and contamination in prokaryotic genomes resulting from misbinning of contigs of unrelaated lineages 

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

In your own words (2 sentences max), explain what is a chimeric bin.


# Manual Bin refinement 
### Anvirefine
- Lets us work on bins manually

In Terminal 

```
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

anvi-refine -c ../5_anvio_profiles/contigs.db -C METABAT -p ../6_anvimerge/PROFILE.db --bin-id METABAT__6

anvi-refine -c ../5_anvio_profiles/contigs.db -C METABAT -p ../6_anvimerge/PROFILE.db --bin-id METABAT__10

```
METABAT 6 looks good, no bins deleted

METABAT 10 overall coverage too low, sample is discarded

## Coverage visualisation: 

```anvi-inspect -p ../6_anvimerge/PROFILE.db -c ../5_anvio_profiles/contigs.db --split-nam ```

```anvi-script-get-coverage-from-bam -b  -C collection-txt -m bin```

## Day 5 

## Taxonomic assignment 

anvi-run-scg-taxonomy associates the single-copy core genes in your contigs-db with taxnomy information

```anvi-run-scg-taxonomy -c ./5_anvio_profiles/contigs.db -T 20 -P 2```

Now you can run anvi-estimate-scg-taxonomy, ‘This program makes quick taxonomy estimates for genomes, metagenomes, or bins stored in your contigs-db using single-copy core genes. 

To estimate abundance of Ribosomal RNAs within your dataset (coverage) use:

```anvi-estimate-scg-taxonomy -c ./5_anvio_profiles/contigs.db -p ../6_anvimerge/PROFILE.db --metagenome-mode --compute-scg-coverages --update-profile-db-with-taxonomy > temp.txt```

ONE final summary to get comprehensive info about your METABAT2 bins:

```anvi-summarize -p ./6_anvimerge/PROFILE.db -c ./5_anvio_profiles/contigs.db --metagenome-mode -o ./SUMMARY_METABAT2 -C METABAT2```




