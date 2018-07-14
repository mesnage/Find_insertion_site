This markdown contains the code used to characterise the insertion site of a lentivirus vector using a whole-genome sequencing approach.


## Alignment to the reference sequence

First of all, the reads are aligned to the 3’LTR sequence using BBmap (Bushnell B. (2015). BBMap. Available at: http://sourceforge.net/projects/bbmap/). Importantly, we use a ‘minratio’ of 0.2 in order to retain the soft-clipping reads which align to only 20% of the target sequence.  


```{bash}

# Align the reads to the insert sequence retaining in soft-clipping reads which aligned to 20% of the target sequence

bbmap.sh minratio=0.2 in1=forward.fq.gz in2=reverse.fq.gz ref=LTR.fa out=mapped.sam
samtools view -S -b mapped.sam > mapped.bam
samtools flagstat mapped.bam
samtools view -b -F 4 mapped.bam > mapped_sorted.bam


```

## Import the BAM file in R to filter the soft clipping reads and produce a table of the results

Reads are extracted from the resulting BAM file using a regular expression searching for the pattern ‘S’ in their cigar strings in R. 

```{r BAM, echo=FALSE}

#Convert BAM file to dataframe using the R package samtools
require(Rsamtools)
bam_file <- scanBam("mapped_sorted.bam")
lst <- lapply(names(bam_file[[1]]), function(elt) {do.call(c, unname(lapply(bam_file, "[[", elt)))})
names(lst) <- names(bam_file[[1]])
bam_df <- do.call("DataFrame", lst)

#use a regular expression to extract soft-clipping reads
bam_df <- bam_df[grep("[S]", bam_df$cigar),]

#save the final results in a table
write.table(bam_df, "soft-clipping_reads.txt", sep ='\t')

```

Ultimately, the soft clipping reads are re-aligned to the 3’LTR sequence using Clustal Omega (https://www.ncbi.nlm.nih.gov/pubmed/25501942) in order to isolate the reads with the longest overhanging sequence. 

This consensus overhanging sequence is then aligned to the human genome using BLAT (https://www.ncbi.nlm.nih.gov/pubmed/11932250) in order to find the location of the insertion site.

Copyright:  (c) Robin Mesnage, King's College London, 2018
Author:     Robin Mesnage
Contact:    robin.mesnage@kcl.ac.uk
