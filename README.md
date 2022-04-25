# ATAC_network
A method to select ATAC peaks which may involve the most cell activities through a correlation network, and check whether those chosen peaks associate with clinical characters.

**1) Select ATAC peak nodes with the most connections in a correlation network**
```
Rscript peak.cor.R TCGA-ATAC_PanCan_Log2Norm_Counts.example.txt > cor.txt
```
This step will take very long time, you can also specify the range of rows needed to be calculated in each session and start several sessions to reduce total time cost. For example:
```
Rscript peak.cor.R TCGA-ATAC_PanCan_Log2Norm_Counts.example.txt 1 100 > cor.1-100.txt
Rscript peak.cor.R TCGA-ATAC_PanCan_Log2Norm_Counts.example.txt 101 200 > cor.101-200.txt
# continue
```
After this step, please select the peaks with top 10% correlations by hand.

**2) Perform principal component analysis on those chosen peaks and check whether some principal component associate with clinical characters, e.g., smoking.**

```
Rscript analysis.R
```
Because analysis.R will load data from MEAN.txt (you could replace this file with data of peaks selected in previous step), please decompress the file MEAN.tar.gz which contains chosen peaks before this step.
