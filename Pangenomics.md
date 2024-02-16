# Our own Pangenome 

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

![My_pangenome_page-0001](https://github.com/MarieRaasch/biol217/assets/157317805/b3dcf94b-7f64-495b-ad37-0c488d191949)


[My_pangenome.pdf](https://github.com/MarieRaasch/biol217/files/14301483/My_pangenome.pdf)


