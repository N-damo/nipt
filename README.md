# TsiOne
## Description  
Three paired chromosomes(chr13,chr18,chr21) were constructed for comparison reference. After GC and mappability correction by LOESS,each bin's read counts was corrected.Then,wilcox rank test was performed to test the significant difference between each pairs.  
## Usage
`python nipt.py -bin 50000 -b *bam -c hg38.chroms.size -m k100.umap.bedgraph -g loess.R -w wilcox.R`  
## Install  
`git clone https://github.com/N-damo/nipt.git`