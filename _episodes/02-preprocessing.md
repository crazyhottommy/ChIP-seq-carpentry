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

### get the raw fastq data

```bash
# go to your home directory
cd ~
# make a new directory
mkdir ChIP-seq

## go insde
cd ChIP-seq

## coopy the fastqs
cp /courses/rai_ChIP-seq/*fq .

## have a look
ls

```


#### quality control

```bash
fastqc IMR90_H3K4me3_chr6.fq
fastqc IMR90_Input_chr6.fq
```

#### Alignment

I will walk through you for the H3K4me3 IP fastq file.

```bash
# use only 1 cpu
bowtie --chunkmbs 320 -m 1 --best -p 1 /courses/bowtie_index/hg19 -q IMR90_H3K4me3_chr6.fq -S > IMR90_H3K4me3_chr6.sam
### reads processed: 247026
# reads with at least one reported alignment: 246922 (99.96%)
# reads that failed to align: 46 (0.02%)
# reads with alignments suppressed due to -m: 58 (0.02%)
# Reported 246922 alignments to 1 output stream(s)

#real    1m30.854s
#user    1m28.444s
#sys     0m1.181s

### save the standard error
bowtie --chunkmbs 320 -m 1 --best -p 1 /courses/bowtie_index/hg19 -q IMR90_H3K4me3_chr6.fq -S > IMR90_H3K4me3_chr6.sam 2> bowtie.log


## check the sam file
less -S IMR90_H3K4me3_chr6.sam

## conver the sam to bam, bam is a binary form of sam
samtools view -Sb -F 4 IMR90_H3K4me3_chr6.sam > IMR90_H3K4me3_chr6.bam

## remove duplicates (there is no duplicates in this example)
## remove duplicates that have exactly the same start and end coordinates. most likely
## due to PCR over-amplification
## -s for single end; -S for paired-end
samtools rmdup -s IMR90_H3K4me3_chr6.bam IMR90_H3K4me3_chr6_rmdup.bam

## sort the bam by coordinates

samtools sort -m 2G -@ 1 -T IMR.tmp -o IMR90_H3K4me3_chr6_rmdup.sorted.bam IMR90_H3K4me3_chr6_rmdup.bam

## index the bam
samtools index IMR90_H3K4me3_chr6_rmdup.sorted.bam

# IMR90_H3K4me3_chr6_rmdup.sorted.bam index will be created.

# check again
ls

## view the alignments
samtools view -h IMR90_H3K4me3.sorted.bam | less -S
```

This step by step process is OK, but it generates too many intermediate files.

For the power users, we use `|` pipe to chain all the step together:

```bash

bowtie --chunkmbs 320 -m 1 --best -p 1 /scratch/genomic_med/apps/annot/indexes/bowtie/hg19 -q IMR90_H3K4me3_chr6.fq -S |  samtools view -Sb -F 4 - | samtools rmdup -s /dev/stdin /dev/stdout |  samtools sort -m 2G -@ 1 -T IMR.tmp -o IMR90_H3K4me3_chr6_rmdup.sorted.bam

samtools index IMR90_H3K4me3_chr6_rmdup.sorted.bam


bowtie --chunkmbs 320 -m 1 --best -p 1 /scratch/genomic_med/apps/annot/indexes/bowtie/hg19 -q IMR90_Input_chr6.fq -S |  samtools view -Sb -F 4 - | samtools rmdup -s /dev/stdin /dev/stdout |  samtools sort -m 2G -@ 1 -T Input.tmp -o IMR90_Input_chr6_rmdup.sorted.bam
```

#### get statistics of the bam file

```bash
samtools flagstat IMR90_H3K4me3_chr6.bam
samtools flagstat IMR90_H3K4me3_chr6_rmdup.sorted.bam
```

#### your turn to align the Input file


#### Peak calling

peak calling without Input control

```bash

macs14 -t IMR90_H3K4me3_chr6_rmdup.sorted.bam -f BAM -g hs --outdir peaks -n IMR90_H3K4me3_no_Input -p 1e-5 --bdg
```

peak calling with Input control

```bash
macs14 -t IMR90_H3K4me3_chr6_rmdup.sorted.bam -c IMR90_Input_chr6_rmdup.sorted.bam -f BAM -g hs --outdir peaks -n IMR90_H3K4me3_with_Input -p 1e-5 --bdg
```

#### bedtools to compare the peak sets

```bash

bedtools intersect -a
```

#### visualize in IGV

peaks and bedgraphs are the two files that you will need to download to your local computer for IGV visualization.

** go to your own local computer**

```bash
rsync -avhP username@101.59:~/ChIP-seq/   .
```
