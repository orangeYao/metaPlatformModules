# metaPlatformModules

## funcPredict.py

### Description:

funcPredict takes in a file C containing a matarix showing functional corelation between a group of genes together with file B with a list of genes with known functions. The resulting output gives a prediction of genes whose functions remain unknown.  

### Usage:
	python functionPrediction/funcPredict.py -c sampleCorr.csv -b sampleB.csv -o sampleOut.csv 
	(default to read from ../../data without flag indicated)

### Arguments: 
	-c: file containing a square matrix where both colum and row names represent gene names, element Aij located at ith row and jth column is set to 1 if there exists known functional relationship between gene i and gene j, otherwise it's set to zero.

	-b: file containing a list of gene names whose function is already known to us

	-o: indicate output name (default as sampleOut.csv)	




## cutTree_r.py

### Description:

cutTree_r takes in a file E containing data about espression regarding to different samples and genes, then output Topological Overlap Matrix and a branch pruning of hierarchical clustering dendrograms. R package "WGCNA" is required.

### Usage: 
	python cut/cutTree_r.py -e sampleEspression.csv -t sampleOutDissTom.csv -c sampleOutColorh.csv
	(default to read from ../../data without flag indicated)

### Arguments:
	-e: file containing espression of genes for different samples, with genes as row names and samples as column names 
	
	-t: indicate output name for diss Topological Overlap Matrix (default as sampleOutDissTom.csv)
	
	-c: indicate output name for dynamic tree cut (default as sampleOutColorh.csv)


## py_svaBatch.py

### Description:

py_svaBatch takes in a file E containing data about espression regarding to different samples and genes, and a file I containging a description including concerned features regarding to different samples. R package "sva" is required. A visualization of such operation by clustering and PCA relying on R package "dendextend" will be accesible soon.

### Usage:
	python batch/py_svaBatch.py -e sampleEspression.csv -i sampleInfo.csv -o sampleOut.csv 

### Arguments:
	-e: file containing espression of genes for different samples, with genes as row names and samples as column names

	-i: file containing information about concerned features regarding to different samples

	-o: indicate output name (default as sampleOut.csv)
