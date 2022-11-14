##
##  PolyGyration.py
##  LatticePoly
##
##  Created by mtortora on 02/08/2020.
##  Copyright © 2020 ENS Lyon. All rights reserved.
import os
import sys

import numpy as np

from scipy.spatial import cKDTree

from vtkReader import vtkReader

from scipy.spatial.distance import pdist, squareform
import time

class Mixing():
	
	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
					
		self.diffRcmFile = os.path.join(self.reader.outputDir, str(time.time())+"cis_trans_duringRepl.res")

		self.Nchain=0
		for t in range(self.reader.nTad):
			if(self.reader.status[t]==-1 or self.reader.status[t]==0):
				self.Nchain+=1


	def Compute(self):
		self.diff=[]
		self.n_mon=[]

		for i in range(self.reader.N):
			self.ProcessFrame(i)
			
			#if (i+1) % 10 == 0:
			#	print("Processed %d out of %d configurations" % (i+1, self.reader.N))
			
				
	def ProcessFrame(self, i):
		data = next(self.reader)
		if(data.nTad<2*self.Nchain and data.nTad>self.Nchain):
			inter=0
			intra=0
			tree1	= cKDTree(data.polyPos[:], boxsize = None)
			pairs = tree1.query_pairs(r = r*0.71)
			for (i,j) in pairs:
				if(data.status[i]==data.status[j] and data.status[i]!=0):
					intra+=1
				if (data.status[i]!=data.status[j] and data.status[i]!=0 and data.status[j]!=0):
					inter+=1
			
			n=data.nTad-self.Nchain
			if(n!=1):
				intra=intra/(n*(n-1))
			inter=inter/n**2
			if(inter!=0):
				self.diff.append(intra/inter)
			else:
				self.diff.append(np.nan)
			self.n_mon.append(data.nTad-self.Nchain)



	def Print(self):
		cumul=[self.n_mon,np.array(self.diff)]
		cumul=np.array(cumul)
		np.savetxt(self.diffRcmFile,cumul.T)
		

		
if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("\033[1;31mUsage is %s outputDir initFrame r\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	r=int(sys.argv[3])

	mix = Mixing(outputDir, initFrame=initFrame)

	mix.Compute()
	mix.Print()
