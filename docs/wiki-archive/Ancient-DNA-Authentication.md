Ancient DNA (aDNA) studies provide invaluable insights into past organisms, environments, and populations. However, aDNA is often highly degraded, fragmented, and contaminated with modern DNA, posing significant challenges for reliable analysis. Authentication protocols address these issues by assessing key indicators such as fragment length, damage patterns (e.g., cytosine deamination at fragment ends), and the presence of environmental or contextual DNA. Strict laboratory practices, including dedicated clean rooms and contamination controls, are essential. Together, these measures ensure that the recovered DNA sequences are genuinely ancient and scientifically valid.

## Degree of fragmentation 

Read lengths of all aligned reads can give you an idea of the overall degree of fragmentation of your data. Also, it can be important in later downstream analysis to rerun on bam files mapped to single species in-relation to understand whether you have DNA from modern and ancient or in rare cases both. 

read_length_distribution plot. 

```{bash}

Extract read lengths from the BAM file

samtools view "$bamfile" | awk '{print length($10)}' > "${sample_name}_read_lengths.txt"

Rscript -e "
    library(ggplot2);
    read_lengths <- read.table('${sample_name}_read_lengths.txt', header = FALSE, col.names = c('Length'));
    plot <- ggplot(read_lengths, aes(x = Length)) +
        geom_histogram(binwidth = 1, fill = 'skyblue', color = 'black') +
        labs(title = '${sample_name} Read Length Distribution', x = 'Read Length (bp)', y = 'Frequency') +
        theme_minimal();
    ggsave('${sample_name}_read_length_distribution.pdf', plot = plot, width = 8, height = 6)

```

For the ease of this, we have made a small script that can be found in the script folder, which loops over all bam files
```{bash}
scripts/readlengthplot_from_bam.sh 
```


## Extract damaged reads
Another important set of tools that have been developed for particularly extracting the damaged reads from a set of aligned reads are (PMDtools)[https://github.com/pontussk/PMDtools] and  

First we run (ngsBriggs)[https://github.com/RAHenriksen/ngsBriggs] which uses a fast maximum likelihood method to estimate post-mortem damage, and writes this likelihood to the MD: score in the bam file. This is a tools that was developed for ancient human genomic datasets, and thus takes into account the modern background contamination fraction. here it is set to 30% (-eps). 
```{bash}
ngsbriggs -bam {}.bam -model nb -eps 0.3 -ibam filein.bam -obam fileout.MD.bam -isrecal 1 -nthread 4
```


```{bash}
scripts/split_bam_by_PMD.sh
```
