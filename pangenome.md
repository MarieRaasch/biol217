# Our own Pangenome 

Creating own directory for our own pangenome: 
`mkdir pangenomics2`

## 1. 

```
cd $WORK/pangenomics/pangenomics2/Genomes_pangenomics

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
done```

## Creating pangenome 

### Visualise 

## Create Pangenome 

```
anvi-gen-genomes-storage -e external-genomes.txt -o My-GENOMES.db

anvi-pan-genome -g My-GENOMES.db --project-name My_pangenome --num-threads 4  

```
anvi-display-pan -p My_pangenome/My-GENOMES.db -g My-GENOMES.db

anvi-compute-genome-similarity -e external-genomes.txt -o ANI -p ./My_pangenome/My_pangenome-PAN.db -T 12