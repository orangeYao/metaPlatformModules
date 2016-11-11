import numpy as np
import pandas as pd
import sys, os
import csv
import time
import scipy
from scipy import linalg
from pandas import Series
sys.path.append(os.path.join(os.path.dirname(__file__),'../../../'))
from NetworkAnalysisPackage.database import Raw_Data

from NetworkAnalysisPackage.parameter import path2abundance
from NetworkAnalysisPackage.parameter import path2sampleInfo
from NetworkAnalysisPackage.parameter import numOfTrimmedData
from NetworkAnalysisPackage.preprocessing import filterByNumOfZeros

from numpy.linalg import inv
from NetworkAnalysisPackage.functions3 import conjugateGradient
from NetworkAnalysisPackage.functions3 import diagonalMatrix
from NetworkAnalysisPackage.functions3 import calculateMatrixA 
from NetworkAnalysisPackage.functions3 import nodeLabelBias

from optparse import OptionParser
parser = OptionParser()

def funcPredict(filenameCorr, filenameB):
	if filenameCorr != None:
		corrcoef = pd.read_csv(filenameCorr,index_col=0)
	else:
		data = Raw_Data(path2abundance)
		info = Raw_Data(path2sampleInfo)
		df1 = filterByNumOfZeros(data.df, numOfTrimmedData)
		corrcoef = diagonalMatrix(df1)

	##1. to find A = I - L = I -(D-W) = I - D + W
	dimension = corrcoef.shape[0]
	np_a = calculateMatrixA(corrcoef, dimension)

	##2. to find b from input list
	if filenameB != None:
		with open(filenameB, 'r') as f:
			gene_list = [int(line.rstrip('\n')) for line in f]
	else:
		gene_list = corrcoef.index.tolist()[20:25]  #simulate input list temporarily

	nega_gene_list = []
	label_bias_b = nodeLabelBias(gene_list, nega_gene_list, dimension, corrcoef.index)

	##3. to find inv(A)*b with diagonalMatrix algorithm
	if np.linalg.det(np_a) != 0:
		label_prop_score2 = np.linalg.solve(np_a, label_bias_b)
		df_out = pd.DataFrame(label_prop_score2,index = corrcoef.index.tolist())
		return df_out

	else:
		print "pseudo inverse of A is used instead"
		label_prop_score2 =  np.linalg.pinv(np_a).dot(label_bias_b)
		df_out = pd.DataFrame(label_prop_score2,index = corrcoef.index.tolist())
		return df_out

if __name__ == "__main__":
	parser.add_option("-c", "--fileCorr", dest="filenameCorr", help="open correlation file", metavar="FILE")
	parser.add_option("-b", "--fileB", dest="filenameB", help="open functional file", metavar="FILE")
	parser.add_option("-o", "--fileOut", dest="filenameOut", help="output file name", metavar="FILE", default="sampleOut.csv")
	(options, args) = parser.parse_args()

	df_out = funcPredict(options.filenameCorr, options.filenameB)
	print df_out
	df_out.to_csv(options.filenameOut)


