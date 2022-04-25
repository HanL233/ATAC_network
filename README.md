# ATAC_network
A method to select ATAC peaks which may involve the most cell activities through a correlation network, and check whether those chosen peaks associate with clinical characters.

**1) Select ATAC peak nodes with the most connections in a correlation network**
```
Rscript peak.cor.R TCGA-ATAC_PanCan_Log2Norm_Counts.example.txt > cor.txt
```
This step will take very long time, you can also specify range of rows needed to be calculated in each session and start servral sessions to reduce the time cost. For example:
```
Rscript peak.cor.R TCGA-ATAC_PanCan_Log2Norm_Counts.example.txt 1 100 > cor.1-100.txt
Rscript peak.cor.R TCGA-ATAC_PanCan_Log2Norm_Counts.example.txt 101 200 > cor.101-200.txt
```
After this step, please select the peaks with top correlations by hand.

**2) perform principal component analysis on those chosen peaks and check whether some principal component associate with clinical characters, e.g., smoking.**

```
Rscript analysis.R
```
Because analysis.R will load data from MEAN.txt, please decompress the file MEAN.tar.gz which contains chosen peaks befroe this step.
