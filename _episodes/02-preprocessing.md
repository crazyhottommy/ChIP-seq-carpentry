---
title: "ChIP-seq preprocessing"
teaching: 0
exercises: 0
questions:
- "How to align ChIP-seq fastq files?"
- "How do you analyze the ChIP-seq data"
objectives:
- "Use bowtie1 for mapping reads"
- "Use MACS for calling peaks"
keypoints:
- "bowtie1 is used for short reads < 70bp. ChIP-seq is quite commmon with 36bp reads"
- "Peak calling"
---

#### quality control


```bash
fastqc IMR90_H3K4me3_chr6.fq
fastqc IMR90_Input_chr6.fq
```

#### Alignment

I will walk through you for the H3K4me3 IP fastq file.

```bash
bowtie --chunkmbs 320 -m 1 --best -p 1 /courses/bowtie_index/hg19 -q IMR90_H3K4me3_chr6.fq -S > IMR90_H3K4me3_chr6.sam

### save the standard error
bowtie --chunkmbs 320 -m 1 --best -p 1 /courses/bowtie_index/hg19 -q IMR90_H3K4me3_chr6.fq -S > IMR90_H3K4me3_chr6.sam 2> bowtie.log


## check the sam file
less -S IMR90_H3K4me3_chr6.sam

## remove duplicates
IMR90_H3K4me3_chr6.sam | samtools rmdup -s > IMR90_H3K4me3_chr6_rmdup.sam

## conver the sam to bam
IMR90_H3K4me3_chr6_rmdup.sam | samtools view -Sb -F 4 - > IMR90_H3K4me3_chr6_rmdup.bam

## sort the bam

samtools sort -m 2G -@ 1 -T IMR.tmp -o IMR90_H3K4me3_chr6_rmdup.sorted.bam IMR90_H3K4me3_chr6_rmdup.bam

## index the bam
samtools index IMR90_H3K4me3_chr6_rmdup.sorted.bam

# a bai index will be created.

## view the alignments
samtools view -h IMR90_H3K4me3.sorted.bam | less -S
```

This step by step process is OK, but it generates too many intermediate files.

For the power users, we use `|` pipe to chain all the step together:

```bash

bowtie --chunkmbs 320 -m 1 --best -p 1 /courses/bowtie_index/hg19 -q IMR90_H3K4me3_chr6.fq -S | samtools rmdup -s | samtools view -Sb -F 4 - | samtools sort -m 2G -@ 1 -T IMR.tmp -o IMR90_H3K4me3_chr6_rmdup.sorted.bam

samtools index IMR90_H3K4me3_chr6_rmdup.sorted.bam

```

#### your turn to align the Input file

#### Peak calling

```bash
macs1 --bedgraph
```
#### visualize in IGV

peaks and bedgraphs
