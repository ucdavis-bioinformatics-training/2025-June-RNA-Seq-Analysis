#!/bin/bash

#SBATCH --job-name=star_index # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=120
#SBATCH --mem=40000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production
#SBATCH --reservation=rnaseqworkshop
#SBATCH --account=workshop
#SBATCH --output=slurmout/star-index_%A.out # File to which STDOUT will be written
#SBATCH --error=slurmout/star-index_%A.err # File to which STDERR will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=myemail@email.com

start=`date +%s`
echo $HOSTNAME

outpath="References"
mkdir -p ${outpath}

cd ${outpath}
#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/GRCm39.primary_assembly.genome.fa.gz
#gunzip GRCm39.primary_assembly.genome.fa.gz
#FASTA="../GRCm39.primary_assembly.genome.fa"
FASTA="/share/workshop/mrnaseq_workshop/Data/GRCm39.primary_assembly.genome.fa"

#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/gencode.vM37.basic.annotation.gtf.gz
#gunzip gencode.vM37.basic.annotation.gtf.gz
#GTF="../gencode.vM37.basic.annotation.gtf"
GTF="/share/workshop/mrnaseq_workshop/Data/gencode.vM37.basic.annotation.gtf"

mkdir star.overlap100.gencode.M37
cd star.overlap100.gencode.M37

module load star

call="STAR
    --runThreadN ${SLURM_NTASKS} \
    --runMode genomeGenerate \
    --genomeDir . \
    --genomeFastaFiles ${FASTA} \
    --sjdbGTFfile ${GTF} \
    --sjdbOverhang 100"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
