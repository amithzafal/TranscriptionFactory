# -*- coding: utf-8 -*-
##
# MonomerDist.py
# LatticePoly
##
# Based on mtortora's  code.
# Copyright © 2021 ENS Lyon. All rights reserved.
##

import os
import sys
import pandas as pd
import numpy as np
import cooler
from scipy.spatial import cKDTree
from cooler.create import ArrayLoader
from vtkReader_multi import vtkReader

from scipy.spatial.distance import pdist, squareform


class MonomerDmap():
	def __init__(self, outputDir, initFrame):
		self.readers=[]
		for i in range(0,16):
			self.readers.append(vtkReader(outputDir, i,initFrame,readLiq=False, readPoly=True))
		self.contactFile = os.path.join(outputDir,"r_"+str(r)+"_"+str(initFrame)+ "_"+str(interval)+"full_genome.cool")


		#self.Nchain=0
		#for t in range(self.reader.nTad):
		#	if(self.reader.status[t]==-1 or self.reader.status[t]==0):
		#		self.Nchain+=1
		#self.timepoint=initFrame #find frame where replication reach the desired percentage
		#if(self.reader.nTad>self.Nchain):
		#	print("Chromosome already replicating ")
		#self.timepoint=init_time #find frame where replication starts
		
				

		#print(self.timepoint)
		#self.reader = vtkReader(outputDir, self.timepoint,readLiq=False, readPoly=True)
		#restarted vtk reader from middle frame of desired percentage
		#compute the hic for the minutes
		self.Compute(interval)

		self.Print()
			






		
	#NB here is not the finalFrame but the number of iterations
	def Compute(self,finalFrame):
		#self.polyAniso = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)
		n_bins=0
		for reader in self.readers:
			n_bins+=len(reader.polyPos)
		self.contactProb = np.zeros((n_bins, n_bins), dtype=np.float32)

		

		for i in range(0, finalFrame):
			self.ProcessFrame(i)
			if (i+1) % 10 == 0:
				print("Processed %d out of %d configurations" %
					  (i+1, finalFrame))


#self.contactProb=np.rint(self.contactProb/finalFrame)

	def ProcessFrame(self, i):
		full_positions=[]
		for reader in self.readers:
			data = next(reader)
			full_positions.append(data.polyPos)

				 
		tree1	= cKDTree(np.concatenate(full_positions), boxsize = None)
		pairs = tree1.query_pairs(r = r*0.71) # NN distance FCC lattice 1/np.sqrt(2) = 0.71
		for (i,j) in pairs:
			self.contactProb[i,j] = self.contactProb[i,j] + 1
			self.contactProb[j,i] = self.contactProb[j,i] + 1
			
			

					



	def Print(self):
		Nchains=[]
		for reader in self.readers:
			Nchains.append(1000*len(reader.polyPos))
		chrom_names=["chr" + str(i+1) for i in range(len(Nchains))]
		chromsizes=pd.Series(Nchains)
		chromsizes=chromsizes.rename(lambda x: chrom_names[x])
		chromsizes=chromsizes.astype('int64')	
		bins = cooler.binnify(chromsizes, 1000)
		pixels = ArrayLoader(bins, self.contactProb, chunksize=10000000)
		cooler.create_cooler(self.contactFile,bins,pixels)
		
		
		print("\033[1;32mPrinted avg.contact probability to '%s'\033[0m" %self.contactFile)


if __name__ == "__main__":
	if len(sys.argv) != 5:
		print("\033[1;31mUsage is %s outputDir initFrame r interval \033[0m" % sys.argv[0])
		sys.exit()
	
	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	r=float(sys.argv[3])
	interval=int(sys.argv[4])
	
	#init_time=int(sys.argv[4])


monom = MonomerDmap(outputDir, initFrame=initFrame)

