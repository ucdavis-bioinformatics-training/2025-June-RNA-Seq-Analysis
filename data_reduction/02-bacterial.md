# Bacterial RNASeq data reduction workflow

## The challenge/difference

The major difference between RNASeq in prokaryotes and in eukaryotes is the gene structure. In generally, there are no introns in prokaryotic genes, therefore no splicing.


<p align = "center">
<img src="https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2025-June-RNA-Seq-Analysis/master/data_reduction/alignment_mm_figures/bacterialrnaseq_figures1.png" alt="prokaryotic gene structure" width="700px"/>
</p>

<p align = "center" style="font-family:Times;font-size:12px;">
https://study.com/skill/practice/contrasting-the-regulation-of-gene-expression-in-prokaryotic-eukaryotic-organisms-questions.html
</p>

This intron-less structure in prokaryotic genes dictates how the RNASeq data should be processed, especially at the alignment stage.

  1. One option is to use an alignment tool/aligner that allows for continuous mapping of the sequencing reads to the reference genome: [__Bowtie2__](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml). In the default mode (end-to-end mode), Bowtie2 searches for alignments that matches the sequencing reads from end to end. In another words, the alignments do not allow any trimming or clipping of bases from the sequencing reads. After the alignment, the quantification of gene expression can be done using [__featureCounts__](https://academic.oup.com/bioinformatics/article/30/7/923/232889), by supplying the gene annotation information.

  1. The other option is to use an aligner that aligns/pseudo-aligns the sequencing reads to the reference transcriptome: such as salmon.


1. An example of processing script for RNASeq data in prokaryotes. This assumes the sequencing data has gone through quality control: adapter trimming, quality trimming, length filtering,...

    <pre class="prettyprint"><code class="language-py" style="background-color:333333">#!/bin/bash


    #SBATCH --nodes=1
    #SBATCH --ntasks=8
    #SBATCH --time=360
    #SBATCH --mem=16000 # Memory pool for all cores (see also --mem-per-cpu)
    #SBATCH --partition=production
    #SBATCH --array=1-54
    
    start=`date +%s`
    echo $HOSTNAME
    echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
    
    sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`
    annotation="References/gencode.vM26.annotation.gtf"
    outpath='02-Bowtie2'
    
    echo "SAMPLE: ${sample}"
    
    module load bowtie2

    bowtie2 --very-sensitive -p 12 -x $refP/reference -1 R1.fastq -2 R2.fastq |samtools view -bh -@ 2 -m 5G -o $outpath/${sample}.bam -
    
    samtools index -b ${outpath}/${sample}.bam
    
    module load subread
    
    featureCounts -T 4 --verbose -s 1 \
                  -a ${annotation} \
                  -t gene \
                  -o ${outpath}//${sample}_featurecounts.txt \
                  ${outpath}/${sample}.bam \
                  > ${outpath}/${sample}_featurecounts.stdout \
                  2> ${outpath}/${sample}_featurecounts.stderr
    
    
    end=`date +%s`
    runtime=$((end-start))
    echo $runtime
    </code></pre>

