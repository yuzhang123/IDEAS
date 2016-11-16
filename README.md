—————————————————————————————————————————————————

	IDEAS v1.10 for Linux, by Yu Zhang	
—————————————————————————————————————————————————

**[11/16/2016]
I fixed a bug that probably has affected the performance of IDEAS to date, not necessarily a big issue but it has to do with memory access. So please use this updated version. Also, please always be alert of how your data is distributed. The -log2 option of IDEAS will by default adding 1 to the signal. So if you've already processed the data in some way such that values smaller than 1 are still meaningful (i.e., not noise), then you need to do -log2 0.1, for instance, to add 0.1 instead of 1 to signal. The number should be the maximum value that you think would corrspond to noise or the signals you don't care about. I added 1 by default for read count data, as 1 read count does not imply signal. 

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


