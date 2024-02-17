# Day 7 Panaroo Pangenome

## Run Panaroo statistics 

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=15:00:00
#SBATCH --job-name=panaroo
#SBATCH --output=panaroo.out
#SBATCH --error=panaroo.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load micromamba/1.4.2
export MAMBA_ROOT_PREFIX=$HOME/.micromamba
eval "$(micromamba shell hook --shell=bash)"
module load micromamba/1.4.2
micromamba activate 08_panaroo
#creata a folder for panaroo
mkdir -p $WORK/pangenomics/01_panaroo
# run panaroo
panaroo -i $WORK/pangenomics/gffs/*.gff -o $WORK/pangenomics/01_panaroo/pangenomics_results --clean-mode strict -t 12
micromamba deactivate
module purge
jobinfo
```

summary statistics.txt

Core genes	(99% <= strains <= 100%)	0
Soft core genes	(95% <= strains < 99%)	0
Shell genes	(15% <= strains < 95%)	9849
Cloud genes	(0% <= strains < 15%)	0
Total genes	(0% <= strains <= 100%)	9849

Uploaded the table gene_presence_absence.Rtab into R Studio: 

```
df <- read.delim("gene_presence_absence.Rtab", header = TRUE)


columns_to_sum <- c("X241155E", "X241156E", "X241157E", "X241158E", "X241159E")

# Sum up selected columns
column_sums <- colSums(df[columns_to_sum])
# Calculate Jaccard similarity for each pair of genomes
similarities <- list()

for (i in 1:(ncol(df) - 1)) {
  for (j in (i + 1):ncol(df)) {
    genome1 <- df[, i]
    genome2 <- df[, j]
    intersection <- sum(genome1 & genome2)  # Count the number of genes present in both genomes
    union <- sum(genome1 | genome2)         # Count the number of unique genes in both genomes
    similarity <- intersection / union      # Calculate Jaccard similarity coefficient
    pair <- paste(names(df)[i], names(df)[j], sep = " - ")
    similarities[[pair]] <- similarity
  }
}

# Find the pair with the highest similarity
most_similar_pair <- names(similarities)[which.max(unlist(similarities))]
similarity_score <- similarities[[most_similar_pair]]

print(paste("The most similar pair of genomes is", most_similar_pair, "with a similarity score of", similarity_score))


```
Numer of genes present: 

X241155E: 3534

X241156E: 3935

X241157E: 0 -> No core genes because we have a hit of 0 here !!!

X241158E: 4702

X241159E: 3932

Output of the similarity test in R: "The most similar pair of genomes is X241156E - X241159E with a similarity score of 0.998729674796748"


-> No core genes in out sample: might be artefact 

Visualising the Pangenome in Cytoscape: /work_beegfs/sunam236/pangenomics/01_panaroo/pangenomics_results/final_graph.gml

![final_graph gml](https://github.com/MarieRaasch/biol217/assets/157317805/f599198e-d38d-4772-8bb5-804e940d1412)


# Day 7b Anvio Pangenome - Example data

```
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

# create new folder
mkdir $WORK/pangenomics/02_anvio_pangenomics
```

Downloading data: 
```
curl -L https://ndownloader.figshare.com/files/28965090 -o V_jascida_genomes.tar.gz
tar -zxvf V_jascida_genomes.tar.gz
ls V_jascida_genomes
```

```
cd $WORK/pangenomics/02_anvio_pangenomics/V_jascida_genomes/

ls *fasta | awk 'BEGIN{FS="_"}{print $1}' > genomes.txt

# remove all contigs <2500 nt
for g in `cat genomes.txt`
do
    echo
    echo "Working on $g ..."
    echo
    anvi-script-reformat-fasta ${g}_scaffolds.fasta \
                               --min-len 2500 \
                               --simplify-names \
                               -o ${g}_scaffolds_2.5K.fasta
done

# generate contigs.db
for g in `cat genomes.txt`
do
    echo
    echo "Working on $g ..."
    echo
    anvi-gen-contigs-database -f ${g}_scaffolds_2.5K.fasta \
                              -o V_jascida_${g}.db \
                              --num-threads 4 \
                              -n V_jascida_${g}
done

# annotate contigs.db
for g in *.db
do
    anvi-run-hmms -c $g --num-threads 4
    anvi-run-ncbi-cogs -c $g --num-threads 4
    anvi-scan-trnas -c $g --num-threads 4
    anvi-run-scg-taxonomy -c $g --num-threads 4
done

```

## Visualize contigs.db

```
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

anvi-display-contigs-stats /path/to.your/databases/*db
```
```
srun --reservation=biol217 --pty --mem=16G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --partition=base /bin/bash

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8_biol217
anvi-display-contigs-stats /path/to.your/databases/*db
```
```
ssh -L 8060:localhost:8080 sunam236@caucluster.rz.uni-kiel.de

ssh -L 8080:localhost:8080 n100

```
```
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

anvi-display-contigs-stats $WORK/pangenomics/02_anvio_pangenomics/V_jascida_genomes/*db

```

```
srun --pty --mem=16G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --partition=base /bin/bash

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8
anvi-display-contigs-stats $WORK/pangenomics/02_anvio_pangenomics/V_jascida_genomes/*db

```
## Create external genome files 

```
anvi-script-gen-genomes-file --input-dir $WORK/pangenomics/02_anvio_pangenomics/V_jascida_genomes/ -o external-genomes.txt
                             
```

## Investigate contamination 

```
cd V_jascida_genomes

anvi-estimate-genome-completeness -e external-genomes.txt > contamination_completeness.txt

```

## Visualize contigs for refinement 

```
anvi-profile -c V_jascida_52.db --sample-name V_jascida_52 --output-dir V_jascida_52 --blank

```

```
srun --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --partition=base /bin/bash

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

anvi-interactive -c V_jascida_52.db -p V_jascida_52/PROFILE.db
                 
```

```
ssh -L 8060:localhost:8080 sunam226@caucluster.rz.uni-kiel.de

ssh -L 8080:localhost:8080 n100
```

```
srun --pty --mem=10G --nodes=1 --tasks-per-node=1 --cpus-per-task=1 --partition=base /bin/bash

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8_biol217

anvi-interactive -c V_jascida_52.db -p V_jascida_52/PROFILE.db
```
## Splitting the genome in our good bins 

```
anvi-split -p V_jascida_52/PROFILE.db -c V_jascida_52.db -C default -o V_jascida_52_SPLIT

# V_jascida_52_SPLIT/V_jascida_52_CLEAN/CONTIGS.db

sed 's/V_jascida_52.db/V_jascida_52_SPLIT\/V_jascida_52_CLEAN\/CONTIGS.db/g' external-genomes.txt > external-genomes-final.txt

```

## Compute Pangenome

```
anvi-estimate-genome-completeness -e external-genomes.txt
anvi-estimate-genome-completeness -e external-genomes-final.txt
```
-> V_jascida_52 hat nun auch eine redundancy von 5.63% 

```
anvi-gen-genomes-storage -e external-genomes-final.txt -o V_jascida-GENOMES.db

anvi-pan-genome -g V_jascida-GENOMES.db --project-name V_jascida --num-threads 4                         
 ```

## Display Pangenome

```
anvi-display-pan -p V_jascida/V_jascida-PAN.db -g V_jascida-GENOMES.db
```

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=5:00:00
#SBATCH --job-name=anvio_pangenomics
#SBATCH --output=anvio_pangenomics.out
#SBATCH --error=anvio_pangenomics.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd $WORK/pangenomics/pangenomics2/Genomes_pangenomics

# ls *fasta > genomes.txt

# remove all contigs <2500 nt
for g in `cat genomes.txt`
do
    echo
    echo "Working on $g ..."
    echo
    anvi-script-reformat-fasta ${g}.fasta \
                               --min-len 2500 \
                               --simplify-names \
                               -o ${g}_scaffolds_2.5K.fasta
done

# generate contigs.db
for g in `cat genomes.txt`
do
    echo
    echo "Working on $g ..."
    echo
    anvi-gen-contigs-database -f ${g}_scaffolds_2.5K.fasta \
                              -o ${g}.db \
                              --num-threads 4 \
                              -n ${g}
done

# annotate contigs.db
for g in *.db
do
    anvi-run-hmms -c $g --num-threads 4
    anvi-run-ncbi-cogs -c $g --num-threads 4
    anvi-scan-trnas -c $g --num-threads 4
    anvi-run-scg-taxonomy -c $g --num-threads 4
done
```

### Display

```module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

anvi-display-contigs-stats $WORK/pangenomics/pangenomics2/Genomes_pangenomics/*db```

### External genome files 

```anvi-script-gen-genomes-file --input-dir $WORK/pangenomics/pangenomics2/Genomes_pangenomics/ -o external-genomes.txt
```

#### Investigate contamination 

```cd Genomes_pangenomics
anvi-estimate-genome-completeness -e external-genomes.txt

```

### Visualise 

## Create Pangenome 

```
anvi-gen-genomes-storage -e external-genomes.txt -o My-GENOMES.db

anvi-pan-genome -g My-GENOMES.db --project-name My_pangenome --num-threads 4  
```


# Anvio: Our own Pangenome 

Creating own directory for our own pangenome: 
`mkdir pangenomics2`

## Genomes downloaded from the NCBI database: 
Bacterial Taxon	+ NCBI Reference Sequence

Clostridium botulinum	NC_009495.1

Clostridium perfringens	NZ_CP075979.1

Chlostridium tetani	NC_004557.1

Pseudomonas putida	NC_021505.1

Vibrio cholerae	NZ_CP043554.1

Vibrio harveyi	NZ_CP125875.1

Aliivibrio fischeri	NC_011184.1

Vibrio anguillarum	NZ_CP031479.1

Vibrio parahaemolyticus	NC_004603.1

Vibrio coralliilyticus	NZ_CP063051.1


## 1. 

```sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=5:00:00
#SBATCH --job-name=my_pangenome_refined
#SBATCH --output=my_pangenome_refined.out
#SBATCH --error=my_pangenome_refined.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

cd $WORK/pangenomics/pangenomics2/Genomes_pangenomics

ls *fasta > genomes.txt

# remove all contigs <2500 nt
for g in `cat genomes.txt`
do
    echo
    echo "Working on $g ..."
    echo
    anvi-script-reformat-fasta ${g}.fasta \
                               --min-len 2500 \
                               --simplify-names \
                               -o ${g}_scaffolds_2.5K.fasta
done

# generate contigs.db
for g in `cat genomes.txt`
do
    echo
    echo "Working on $g ..."
    echo
    anvi-gen-contigs-database -f ${g}_scaffolds_2.5K.fasta \
                              -o .db \
                              --num-threads 4 \
                              -n ${g}
done

# annotate contigs.db
for g in *.db
do
    anvi-run-hmms -c $g --num-threads 4
    anvi-run-ncbi-cogs -c $g --num-threads 4
    anvi-scan-trnas -c $g --num-threads 4
    anvi-run-scg-taxonomy -c $g --num-threads 4
done```

## Creating pangenome 

### Visualise 

## Create Pangenome 
```


```sh
anvi-gen-genomes-storage -e external-genomes.txt -o My-GENOMES.db

anvi-pan-genome -g My-GENOMES.db --project-name My_pangenome --num-threads 4  

```
```
anvi-display-pan -p My_pangenome-PAN.db -g My-GENOMES.db

anvi-compute-genome-similarity -e external-genomes.txt -o ANI -p ./My_pangenome/My_pangenome-PAN.db -T 12

```
```sh
anvi-gen-genomes-storage -e external-genomes.txt -o My-GENOMES.db

anvi-pan-genome -g My-GENOMES.db --project-name My_pangenome --num-threads 4  
```

```sh
cd $WORK/pangenomics/pangenomics2

anvi-get-sequences-for-gene-clusters -p ./My_pangenome/My_pangenome-PAN.db -g ./My_pangenome/My-GENOMES.db --min-num-genomes-gene-cluster-occurs 9 --max-num-genes-from-each-genome 1 --concatenate-gene-clusters --output-file ./My_pangenome/My_new_pangenome-SCGs.fa

trimal -in ./My_pangenome/My_new_pangenome-SCGs.fa -out ./My_pangenome/My_new_pangenome-SCGs-trimmed.fa -gt 0.5 

iqtree -s ./My_pangenome/My_new_pangenome-SCGs-trimmed.fa -m WAG -bb 1000 -nt 8

echo -e "item_name\tdata_type\tdata_value" > My_pangenome/My_new_pangenome-phylogenomic-layer-order.txt

# add the newick tree as an order
echo -e "SCGs_Bayesian_Tree\tnewick\t`cat My_pangenome/My_new_pangenome-SCGs-trimmed.fa.treefile`" >> My_pangenome/My_new_pangenome-phylogenomic-layer-order.txt

# import the layers order file
anvi-import-misc-data -p ./My_pangenome/My_pangenome-PAN.db -t My_pangenome/My_new_pangenome-phylogenomic-layer-order.txt
```
![Bildschirmfoto_16-2-2024_162728_](https://github.com/MarieRaasch/biol217/assets/157317805/f232d243-003c-41c1-9bed-2f0a9d53e908)

PDF for better Quality:

[My_pangenome.pdf](https://github.com/MarieRaasch/biol217/files/14312891/My_pangenome.pdf)



