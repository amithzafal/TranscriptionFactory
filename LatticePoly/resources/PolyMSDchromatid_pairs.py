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
from vtkReader import vtkReader
import time



class PolyMSD():

	def __init__(self, outputDir, initFrame):
		self.Nchain=0
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		for t in range(self.reader.nTad):
				if(self.reader.status[t]==-1 or self.reader.status[t]==0):
					self.Nchain+=1
		print(self.Nchain)
		
		self.msdHetFile = os.path.join(self.reader.outputDir,str(time.time())+ "polyHetMSD.res")
		self.msdHomFile = os.path.join(self.reader.outputDir, str(time.time())+"polyHomMSDchromatid1.res")

		if os.path.exists(self.msdHetFile) & os.path.exists(self.msdHomFile):
			print("Files '%s' and '%s' already exist - aborting" % (self.msdHetFile, self.msdHomFile))
			sys.exit()


	def Compute(self):
		vMem = psutil.virtual_memory()
		sizeTot = self.reader.N * self.reader.polyPos.nbytes
		
		if sizeTot < vMem.available:
			self.cumulDistHet = 0
			self.cumulDistHom = 0
		
			posHist = self.ReadHist()
			np.save("posHist",posHist)

			for idxTad in range(self.Nchain):
				if self.reader.polyType[idxTad] == 1:
					self.cumulDistHet += msdFFT(posHist[:, idxTad])
				else:
					pos1=posHist[:, idxTad]
					pos2=posHist[:, self.reader.SisterID[idxTad]]
					diff=pos1-pos2
					self.cumulDistHom += msdFFT(dif)
							
				if (idxTad+1) % 1000 == 0:
					print("Processed %d out of %d TADs" % (idxTad+1, self.reader.nTad))
					
		else:
			print("Memory overflow - reduce chosen number of frames")
			sys.exit()


	def ComputeTad(self, idxTad):
		tadPosHist = np.zeros((self.reader.N, 3), dtype=np.float32)

		for i in range(self.reader.N):
			data = next(self.reader)
			tadPosHist[i] = data.polyPos[idxTad]
			
		self.distTad = msdFFT(tadPosHist)
				

	def ReadHist(self):
		posHist = np.zeros((self.reader.N, self.Nchain, 3), dtype=np.float32)
		
		for i in range(self.reader.N):
			data = next(self.reader)
			posHist[i] = data.polyPos[:self.Nchain]
			
		return posHist
	
	
	def Print(self):
		if self.reader.nHet > 0:
			msdHet = self.cumulDistHet /  self.reader.nHet
			np.savetxt(self.msdHetFile, msdHet)
			
			print("\033[1;32mPrinted heterochromatic MSDs to '%s'\033[0m" % self.msdHetFile)
			
		if self.reader.nEuc > 0:
			msdHom = self.cumulDistHom / self.reader.nEuc
			np.savetxt(self.msdHomFile, msdHom)
			
			print("\033[1;32mPrinted euchromatic MSDs to '%s'\033[0m" % self.msdHomFile)

	
	def PrintTad(self, idxTad):
		msdFile = self.reader.outputDir + "/"
		np.savetxt(msdFile, self.distTad)
		
		print("\033[1;32mPrinted TAD MSD to '%s'\033[0m" % msdFile)
	
	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	msd = PolyMSD(outputDir, initFrame=initFrame)



	if len(sys.argv) == 3:
		msd.Compute()
		msd.Print()
		
	elif len(sys.argv) == 4:
		idxTad = int(sys.argv[3])
	
		msd.ComputeTad(idxTad)
		msd.PrintTad(idxTad)
