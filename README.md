—————————————————————————————————————————————————

	IDEAS v1.20 for Linux, by Yu Zhang	
—————————————————————————————————————————————————

**[06/29/2017]
A major update. The new version 1.2 implements a pipeline that identifies reproducible epigenetic states, and thus greatly improves the prediction robustness. Also, the new version streamlines the analyses of data preprocessing, segmentation, and result visualization together, or separately if desired, within one program. R is required for the new pipelines. Finally, the program can be run in parallel using multiple threads.

----------------------------------------------------------------------------------------------------------------------------------------

**[11/16/2016]
I fixed a bug that probably has affected the performance of IDEAS to date, not necessarily a big issue but it has to do with memory access. So please use this updated version. Also, please be alert of how your data is distributed. The -log2 option of IDEAS will by default add 1 to the signal. If your data signal is mean read count per window, or if you've processed the data in some way such that values smaller than 1 are meaningful (i.e., not noise), then you need to use -log2 0.1, for instance, to add 0.1 instead of 1 to the data when taking log. This number should be the value that you think would corrspond to noise (or the signals below which you won't care about). I added 1 by default for maximum read count per window, as 1 read count does not imply signal. 

**[11/16/2016]
I added a R script, runideaspipe.R, to run ideas in a hopefully more stable way.The motivation is that each time IDEAS can produce slightly different results, especially different number of states, due to different starting values of the model. This R script will automatically start from different start values and combine results from different runs to come up with a concensus, followed by retraining of IDEAS. As a result, the output will be more stable and the states will be more reproducible. It however does not mean that you'll get identical result each time, as the problem that IDEAS solves is not a convex problem. See README_runideaspipe for details.

**[11/16/2016]
I added an option "-hp" for running IDEAS to improve continuity of states along chromosome. This will make the annotation smoothier, but at risk of loosing precision in details. This option is still under testing.

**[11/16/2016]
I added an option '-inv' for running IDEAS in an interval of the genome.

**[08/15/2016] now allow the user to provide BAM or BigWig files as input, or directly download those data from public repositories and automatically process to fed into IDEAS.  

This program is designed to segment genomes in multiple cell types simultaneously so to better identify functional elements and epigenomic variation/conservation patterns, both globally and locally, across all cell types.

See included package for details and usage.

REFERENCES

	Yu Zhang, Feng Yue, Ross C. Hardison. Jointly characterizing epigenetic dynamics across multiple human cell types. 
	Nucleic Acids Res, 44(14):6721-31.

__________________________________________________________
Questions and comments?

please contact me at yzz2 at psu.edu , Thanks,
__________________________________________________________


