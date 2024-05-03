##
##  PolyMSD.py
##  LatticePoly
##
##  Created by mtortora on 15/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys
import psutil

import numpy as np

from utils import msdFFT
from vtkReader_multi import vtkReader
import time

class PolyMSD():

	def __init__(self, outputDir,chrom, initFrame):
		self.reader = vtkReader(outputDir, chrom, initFrame, readLiq=False, readPoly=True)
		
		self.msdHomFile = os.path.join(self.reader.outputDir, str(chrom)+"MSD_after_repl.res")

		
		self.Nchain=1000
		self.timepoint=0 #find frame where replication reach the desired percentage
		
			
			

	def Compute(self):
			self.cumulDistHet = 0
			self.cumulDistHom = 0
			posHist = self.ReadHist()
			for idxTad in range(self.Nchain):
				self.cumulDistHom += msdFFT(posHist[:, idxTad])
					


	
	def ReadHist(self):
		posHist = np.zeros((self.reader.N-self.timepoint, self.Nchain, 3), dtype=np.float32)
				
		for i in range(self.reader.N-self.timepoint):
			data = next(self.reader)
			posHist[i] = data.polyPos
		print("poshist_completed")
		return posHist
	
	
	def Print(self):
		msdHom = self.cumulDistHom / self.Nchain
		np.savetxt(self.msdHomFile, msdHom)
		print("print")			

	

	
	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	chrom=int(sys.argv[3])	
	msd = PolyMSD(outputDir, chrom,initFrame=initFrame)



	if len(sys.argv) == 4:
		msd.Compute()
		msd.Print()
		
