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

# anvi-gen-genomes-storage -e external-genomes.txt -o My-GENOMES.db

# anvi-pan-genome -g My-GENOMES.db --project-name My_pangenome --num-threads 4  

cd $WORK/pangenomics/pangenomics2

anvi-get-sequences-for-gene-clusters -p ./My_pangenome/My_pangenome-PAN.db -g ./My_pangenome/My-GENOMES.db --min-num-genomes-gene-cluster-occurs 9 --max-num-genes-from-each-genome 1 --concatenate-gene-clusters --output-file ./My_pangenome/My_new_pangenome-SCGs.fa

trimal -in ./My_pangenome/My_new_pangenome-SCGs.fa -out ./My_pangenome/My_new_pangenome-SCGs-trimmed.fa -gt 0.5 

iqtree -s ./My_pangenome/My_new_pangenome-SCGs-trimmed.fa -m WAG -bb 1000 -nt 8

echo -e "item_name\tdata_type\tdata_value" > My_pangenome/My_new_pangenome-phylogenomic-layer-order.txt

# add the newick tree as an order
echo -e "SCGs_Bayesian_Tree\tnewick\t`cat My_pangenome/My_new_pangenome-SCGs-trimmed.fa.treefile`" >> My_pangenome/My_new_pangenome-phylogenomic-layer-order.txt

# import the layers order file
anvi-import-misc-data -p ./My_pangenome/My_pangenome-PAN.db -t My_pangenome/My_new_pangenome-phylogenomic-layer-order.txt

