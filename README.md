# IDEAS
1. Installation: 
	Just unzip and untar to a folder.


2. Run:
	./ideas inputfile [options]
	
	Options:
	-log2	take log2(x+1) transformation of the input data, recommented for read count data to reduce skewness.
	-minerr [number]	specify the minimum stardard deviation for each state, default 0.01, but we recommend a larger number, say 0.5 or 1. 
				this option helps to remove artifical states created by tightly concentrated data, such as 0 in read count data.
	-sample [burnin mcmc]	specify the number of burnin and sampling steps. By default it is 50 50.
	
	[this program currently outputs "tmp.out" file, you can just ignore it, unless you are running the job in background and it runs really slow and you want to know which iteration it's at]


3. Input Format:
	1st line: column names
	remaining lines: position information and data

	An example for 3 marks in 2 cell lines:
	chr pos cell1.mark1 cell1.mark2 cell1.mark3 cell2.mark1 cell2.mark2 cell2.mark3
	1 1000200 1 5 0 4 0 12
	1 1000400 4 17 1 1 1 0
	1 1000600 10 31 7 8 2 0
	
	In this example, the first two columns are genome positions, the next 6 columns are read count data.
	The column names indicate that cell 1 data are provided first, and then cell 2 data. 
		"cell1" etc can be replaced by actual cell type names without space.
		"mark1" etc can also be replaced by actual mark names without space.
		both cell type and mark names need to be provided in the 1st line, and they need to be connected by ".".
	The toy.data has idential column names for cell-mark comibnations, which indicate replicate data.
	If certain mark is missing in a cell type, it is fine, such as:
		chr pos cell1.mark1 cell1.mark3 cell2.mark1 cell2.mark2
		1 1000200 1 0 4 0 
		1 1000400 4 1 1 1 
		1 1000600 10 7 8 2 
	In this example, mark2 in cell1 and mark3 in cell2 are missing.
	However, in such case, the program wont know the order of marks as you want them to be (i.e., mark1 mark2 mark3). 
	So the order of marks taken by the program will be determined by their order of appearance, and in this example, mark1, mark3, mark2.
	This matters when read the output from the program, such as the mean vector for each epigenetic state.


4. Output Files:
	*.g file: 	epigenetic states and position classes
			first 3 columns are index, chr, position
			the next N columns are epigenetic states, where N=number of replicates of cell types. For toy.data, N=10 because each of 5 cell types has 2 replicates.
			the last column is the position class label.

	*.para:		frequency, mean and variance parameters for each epigenetic states. 
			read by readPara() function in plot.R

	*.pop file: 	local cell type clustering result, one row for each replicate of cell type. 
			shown by showPop() function in plot.R
	
	*.popr file:	probability of cell type cluster boundary, not very useful.


5. plot.R:
	This is a R script that summarizes results.
	In R environment, load this file by source("plot.R").
	
	get parameters for epigenetic states: readPara(), 
		example: 	readPara("toy.data",1,3)
				first parameter is input file name, 
				second parameter is always 1,
				third parameter is number of marks. if there are missing marks, this number should be the number of marks as if there are no missing marks.
			
				This function outputs an object with 3 components:
				$p: proportion of epigenetic states in the entire data
				$m: mean vector for all epigenetic states, every k values correspond to one state, where k is the number of marks you input
				$v: covariance matrices for all epigenetic states, every k rows correpond to one state

	show cell type clustering and epigenetic states: showPop()
		example:	showPop("toy.data", indlist=1:10, inv=1:1000, repn=2)
				first parameter is input file name,
				second parameter list which (replicate of) cell type should be shown, if not specified, all will be shown,
				third parameter list which positions should be shown, if not specified, all will be shown,
				fouth parameter specifies number of replicates per cell type (only works if the number of replicates for each cell is the same). not needed if no replicates.

				This function will show a figure with three parts, top one is the cell type clustering figure using results in *.pop file. 
				middle one is the epigenetic state using results in *.g file.
				bottom one is the position class using results in *.g file.


6. Questions and comments?
	please contact me at yzz2 at psu.edu, Thanks,
