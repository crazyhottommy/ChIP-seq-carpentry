---
title: "Introduction: ChIP-seq assay"
teaching: 0
exercises: 0
questions:
- "what is ChIP-seq?"
- "How do you analyze the ChIP-seq data"
objectives:
- "Understand what the ChIP-seq assay is measuring"
- "Undertand steps to process raw ChIP-seq fastqs to final peaks for visualization"
- "know what the data set is"
keypoints:
- "ChIP-seq is sequencing specific regions of the genome"
- "ChIP-seq peaks can give insights of where the protein of interest is enriched in the genome"
- "target 20 million reads for ChIP-seq experiment"
---

### What is ChIP-seq?
Chromatin immunoprecipitation (ChIP) followed by high-throughput DNA sequencing (ChIP-seq) is a technique to map genome-wide transcription
factor binding sites and histone-modification enriched regions.

Briefly, DNA bounding proteins and DNA (Chromatin) are cross-linked by formaldehyde and the chromatin is sheared by sonication into small fragments (typically 200 ~ 600 bp). The protein-DNA complex is immnuoprecipitated by an antibody specific to the protein of interest. Then the DNA is purified and made to a library for sequencing. After aligning the sequencing reads to a reference genome, the genomic regions with many reads enriched are where the protein of interest bounds. ChIP-seq is a critical method for dissecting the regulatory mechanisms of gene expression.

### What are the processing steps for ChIP-seq data?
The general processing steps are as following:

1. Quality control of your `fastq` read files using tools such as [FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Read our previous section: [Visualizing sequencing data quality](https://read.biostarhandbook.com/data/fastq-quality-visualization.html).

2. Aligning fastq reads to reference genome. The most popular read aligners are `bwa` and `bowtie`. Read section [Short Read Aligners](https://read.biostarhandbook.com/align/short-read-aligners.html). `bowtie` has two versions: `bowtie1` and `bowtie2`. `bowtie2` is better for reads length greater than 50bp. It is still very common to use 36bp single-end sequencing library for ChIP-seq, so `bowtie1` is preferred for ChIP-seq with short reads.

3. Calling peaks from the alignment bam files.

4. Visualizing the resulting peak files and raw signal files.

### What are IgG control and input control for ChIP-seq?
Like any other experiments, a control is needed for ChIP-seq. There are two kinds of controls for ChIP-seq: IgG control and input control. IgG control is DNA resulting from a "mock" ChIP with Immunoglobulin G (IgG) antibody, which binds to non-nuclear antigen; input control is DNA purified from cells that are cross-linked, fragmented, but without adding any antibody for enrichment.

One problem with IgG control is that if too little DNA is recovered after immunoprecipitation(IP), sequencing library will be of low complexity and binding sites identified using this control could be biased. Input DNA control is ideal in most of the cases. It represents all the chromatin that was available for IP. Read [this biostars post](https://www.biostars.org/p/15817/) for discussion.


### The data set

Download the data, we will use the data from epigenomeRoadmap project.

**The fastq files are already on the server**: `/course/ChIPseq_lab/fastqs`

you do not have to do the following.

ChIP-Seq input data of IMR90 cell line:

https://www.ebi.ac.uk/ena/data/view/SRR037635

H3K4me3 data of the same cell type:

https://www.ebi.ac.uk/ena/data/view/SRR029610

Download the `fastq.gz` file:

IMR90 is a human fibroblast cell line derived from lungs of a 16-week female fetus.

~~~
# H3K4me3
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR029/SRR029610/SRR029610.fastq.gz

# Input
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR037/SRR037635/SRR037635.fastq.gz
~~~
{: .bash}

I extracted the reads only on chromosome 6 to get small size fastqs.
