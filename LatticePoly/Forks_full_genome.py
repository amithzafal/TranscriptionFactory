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
import networkx as nx
from scipy.spatial.distance import pdist, squareform


class forks_full_genome():
	def __init__(self, outputDir, initFrame):
		self.readers=[]
		self.N_chain=[]
		for i in range(0,17):
			self.readers.append(vtkReader(outputDir, i,initFrame,readLiq=False, readPoly=True))
			self.N_chain.append(np.sum([(self.readers[i].status==-1)+(self.readers[i].status==0)]))
		self.ForksFile = os.path.join(outputDir,str(initFrame)+"_forks_full_genome.res")



		self.Compute(self.readers[0].N)

		self.Print()
			






		
	def Compute(self,finalFrame):
		self.Forks_number=[]
		for i in range(0, finalFrame):
			self.ProcessFrame(i)
			if (i+1) % 10 == 0:
				print("Processed %d out of %d configurations" %
					  (i+1, finalFrame))



	def ProcessFrame(self, i):
		forks_for_chrom=[]
		for chrom,reader in enumerate(self.readers):
			data = next(reader)
			fork_mask=np.array((data.fork==1)+(data.fork==-1))
			forks_for_chrom.append(np.sum(fork_mask))
		self.Forks_number.append(np.sum(forks_for_chrom))

				



		
			

					



	def Print(self):
		np.savetxt(self.ForksFile, self.Forks_number)
		


if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s outputDir initFrame  \033[0m" % sys.argv[0])
		sys.exit()
	
	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	#init_time=int(sys.argv[4])


monom = forks_full_genome(outputDir, initFrame=initFrame)

