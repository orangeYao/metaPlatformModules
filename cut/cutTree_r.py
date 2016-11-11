import numpy as np
import pandas as pd
import sys, os
import csv
from pandas import Series
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

sys.path.append(os.path.join(os.path.dirname(__file__),'../../../'))
from NetworkAnalysisPackage.database import Raw_Data

from NetworkAnalysisPackage.parameter import path2abundance
from NetworkAnalysisPackage.parameter import path2sampleInfo
from NetworkAnalysisPackage.parameter import numOfTrimmedData
from NetworkAnalysisPackage.preprocessing import filterByNumOfZeros

from optparse import OptionParser
parser = OptionParser()

WGCNA = importr('WGCNA')
pandas2ri.activate()

def cutTree(filenameEspr):
	if filenameEspr != None:
		df1 = pd.DataFrame.from_csv(filenameEspr)
	else:
		data = Raw_Data(path2abundance)
		df1 = filterByNumOfZeros(data.df, numOfTrimmedData)


	robjects.r(
			'''
			f_dissT <- function(x, beta) {
				allowWGCNAThreads()
				ALLOW_WGCNA_THREADS=48
				options(stringsAsFactors=F)
				datExpr <- t(x)
				ADJ = adjacency(datExpr,power=beta)
				gc()
				dissTOM=TOMdist(ADJ)
				gc()
				return (dissTOM)
			}
			f_cutTree <- function(dissTOM) {
				hierTOM = hclust(as.dist(dissTOM),method="average");
				myheightcutoff =0.55
				mydeepSplit = FALSE # fine structure within module
				myminModuleSize = 2 # modules must have this minimum number of genes
				print (myminModuleSize)
				colorh1=cutreeDynamic(hierTOM, deepSplit=mydeepSplit, cutHeight=myheightcutoff, minClusterSize=myminModuleSize)
				print (myminModuleSize)
				print (hierTOM)
				return (colorh1)
			}
			''')
	r_dissTOM = robjects.r['f_dissT'](df1, 8)
	dissTOM = pandas2ri.ri2py(r_dissTOM)
	r_colorh1 = robjects.r['f_cutTree'](r_dissTOM)
	colorh1 = pandas2ri.ri2py(r_colorh1)
	return dissTOM, colorh1


if __name__ == "__main__":
	parser.add_option("-e", "--fileEs", dest="filenameEspr", help="open correlation file", metavar="FILE")
	parser.add_option("-t", "--fileOutT", dest="filenameOutT", help="output file name", metavar="FILE", default="sampleOutDissTom.csv")
	parser.add_option("-c", "--fileOutC", dest="filenameOutC", help="output file name", metavar="FILE", default="sampleOutColorh.csv")
	(options, args) = parser.parse_args()

	dissTOM, colorh1 = cutTree(options.filenameEspr)	
	np.savetxt(options.filenameOutT, dissTOM, delimiter=",")
	np.savetxt(options.filenameOutC, colorh1, delimiter=",")
