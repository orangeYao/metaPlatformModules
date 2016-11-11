import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

from optparse import OptionParser
parser = OptionParser()
sva = importr('sva')


def svaBatch(filenameEspr, filenameInfo):
	robjects.r(
		'''
		options(stringsAsFactors=F)

		cleanY = function(y, mod, svs) {
			X = cbind(mod, svs)
			Hat = solve(t(X) %*% X) %*% t(X)
			beta = (Hat %*% t(y))
			rm(Hat)
			gc()
			P = ncol(mod)
			return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
		}

		f_sva <- function(dataFileName, infoFileName) {
			setwd(".")
			edata <- read.csv(dataFileName, row.names=1,header=T)
			edata <- as.matrix(edata)
			sample_info <- read.table(file=infoFileName, header=T)

			factorName <- c(colnames(sample_info))
			mod = model.matrix((~as.factor(get(factorName[1]))), data=sample_info)
			mod0 = model.matrix(~1,data=sample_info)

			## to eliminate error, get rid rows with low variance, but seems no use
			getVar <- apply(edata,1,var)
			param <- 1
			data_removeNoVariance <- edata[getVar > param& !is.na(getVar),]

			n.sv = num.sv(edata,mod,method="leek")
			print(n.sv)
			svobj = sva(edata,mod,mod0,n.sv=n.sv-1)
			newEdata <- cleanY(edata,mod,svobj$sv)
			return (newEdata)
		}
		'''
	)

	r_newEdata = robjects.r['f_sva'](filenameEspr, filenameInfo)
	newEdata = pandas2ri.ri2py(r_newEdata)
	return newEdata


def plotClutering():
	dendextend = importr('dendextend')
	robjects.r(
			'''
			'''
			)	


if __name__ == "__main__":
	parser.add_option("-e", "--fileEs", dest="filenameEspr", help="open gene espression file", metavar="FILE")
	parser.add_option("-i", "--fileIn", dest="filenameInfo", help="open sample information file", metavar="FILE")
	parser.add_option("-o", "--fileOut", dest="filenameOut", help="output file name", metavar="FILE", default="sampleOut.csv")
	(options, args) = parser.parse_args()

	newEdata = svaBatch(options.filenameEspr, options.filenameInfo)
	np.savetxt(options.filenameOut, newEdata, delimiter=",")
