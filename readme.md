# Biol 217 Practice Session Day -1

Login Data: ssh -X sunam236@caucluster.rz.uni-kiel.de

What we have learned so far?

1. Basic Linux 
2. Bioinfomrtics basic understandings 
3. Linus comands 

- copy from one folder to another:

Block of code 

```sh
cp source destination 
```

this is the command `cp`

### Task: Practice how to upload images and links 

# Day 2 

## Quality control

Anvio_slurm.txt


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

To run the text `sbatch file.txt`
To check process `squeue -u sunam236`


Kopieren und im Rechner öffnen 

```
scp sunam236@caucluster.rz.uni-kiel.de:/work_beegfs/sunam236/Metagenomics/1_fastqc/*.html .
```


