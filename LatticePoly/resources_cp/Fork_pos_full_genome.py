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


class Fork_pos_full_genome():
	def __init__(self, outputDir, initFrame):
		self.readers=[]
		self.N_chain=[]
		for i in range(0,17):
			self.readers.append(vtkReader(outputDir, i,initFrame,readLiq=False, readPoly=True))
			self.N_chain.append(np.sum([(self.readers[i].status==-1)+(self.readers[i].status==0)]))
		self.Fork_pos = os.path.join(outputDir,str(initFrame)+"_"+str(interval)+"_fork_pos_full_genome.res")




		self.Compute(interval)

		self.Print()
			






		
	def Compute(self,finalFrame):
		self.forks_pos=np.zeros(200)
		for i in range(0, finalFrame):
			self.ProcessFrame(i)
			if (i+1) % 10 == 0:
				print("Processed %d out of %d configurations" %
					  (i+1, finalFrame))



	def ProcessFrame(self, i):
		full_positions=[]
		for chrom,reader in enumerate(self.readers):
			data = next(reader)
			fork_mask=np.array((data.fork==1)+(data.fork==-1))
			if(np.sum(fork_mask)>0):
				full_positions.append(data.polyPos[fork_mask])
		full_positions=np.concatenate(full_positions)
		for pos in full_positions:
			diff=pos-np.array([35,35,70])
			self.forks_pos[round(np.sqrt(np.dot(diff.T,diff))/(2*0.7))]+=1
		


	def Print(self):
		np.savetxt(self.Fork_pos, self.forks_pos)
		

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("\033[1;31mUsage is %s outputDir initFrame r interval \033[0m" % sys.argv[0])
		sys.exit()
	
	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	interval=int(sys.argv[3])
	
	#init_time=int(sys.argv[4])


monom = Fork_pos_full_genome(outputDir, initFrame=initFrame)

