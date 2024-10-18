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
from hdf5Reader import hdf5Reader
import h5py

class PolyMSD():

	def __init__(self, outputDir, fileName, initFrame):
		self.reader = hdf5Reader(outputDir, fileName, initFrame, readLiq=False, readPoly=True)
		self.filePath = os.path.join(outputDir, fileName)
		
		# self.msdHetFile = os.path.join(self.reader.outputDir, "polyHetMSD.res")
		# self.msdHomFile = os.path.join(self.reader.outputDir, "polyHomMSD.res")
		# self.msdPREFile = os.path.join(self.reader.outputDir, "polyPREMSD.res")

		# if os.path.exists(self.msdHetFile) & os.path.exists(self.msdHomFile) & os.path.exists(self.msdPREFile):
		# 	print("Files '%s', '%s' and '%s' already exist - aborting" % (self.msdHetFile, self.msdHomFile, self.msdPREFile))
		# 	sys.exit()


	def Compute(self):
		vMem = psutil.virtual_memory()
		sizeTot = self.reader.N * self.reader.polyPos.nbytes
		
		if sizeTot < vMem.available:
			self.cumulDistHet = 0
			self.cumulDistHom = 0
			self.cumulDistPRE = 0
		
			posHist = self.ReadHist()
			
			for idxTad in range(self.reader.nTad):
				if self.reader.polyPainter[idxTad] == 1:
					self.cumulDistPRE += msdFFT(posHist[:, idxTad])
				elif self.reader.polyType[idxTad] == 1:
					self.cumulDistHet += msdFFT(posHist[:, idxTad])
				else:
					self.cumulDistHom += msdFFT(posHist[:, idxTad])
							
				if (idxTad+1) % 1000 == 0:
					print("Processed %d out of %d TADs" % (idxTad+1, self.reader.nTad))
					
		else:
			print("Memory overflow likely - reduce chosen number of frames")
			sys.exit()


	def ComputeTad(self, idxTad):
		tadPosHist = np.zeros((self.reader.N, 3), dtype=np.float32)

		for i in range(self.reader.N):
			data = next(self.reader)
			tadPosHist[i] = data.polyPos[idxTad]
			
		self.distTad = msdFFT(tadPosHist)
				

	def ReadHist(self):
		posHist = np.zeros((self.reader.N, self.reader.nTad, 3), dtype=np.float32)
		
		for i in range(self.reader.N):
			data = next(self.reader)
			posHist[i] = data.polyPos
			
		return posHist
	
	
	def Print(self):
		self.reader.Close()
		file = h5py.File(self.filePath, 'r+')
		if np.count_nonzero(self.reader.polyPainter == 1) > 0:
			msdPRE = self.cumulDistPRE / np.count_nonzero(self.reader.polyPainter == 1)
			file.create_dataset("process_polyPREMSD", data = msdPRE)

			print("\033[1;32mPrinted PRE MSDs to '%s'\033[0m" % "process_polyPREMSD")

		if self.reader.nHet > 0:
			msdHet = self.cumulDistHet /  self.reader.nHet
			file.create_dataset("process_polyHetMSD", data = msdHet)
			
			print("\033[1;32mPrinted heterochromatic MSDs to '%s'\033[0m" % "process_polyHetMSD")
			
		if self.reader.nEuc > 0:
			msdHom = self.cumulDistHom / self.reader.nEuc
			file.create_dataset("process_polyHomMSD", data = msdHom)
			
			print("\033[1;32mPrinted euchromatic MSDs to '%s'\033[0m" % "process_polyHomMSD")
		
		file.close()

	
	def PrintTad(self, idxTad):
		self.reader.Close()
		file = h5py.File(self.filePath, 'r+')
		msdFile = self.reader.outputDir + "/msdTad%05d.res" % idxTad
		file.create_dataset("process_msdTad%05d" % idxTad, data = self.distTad)
		
		print("\033[1;32mPrinted TAD MSD to '%s'\033[0m" % "process_msdTad%05d" % idxTad)

		file.close()
	
	
if __name__ == "__main__":
	if len(sys.argv) not in [4, 5]:
		print("\033[1;31mUsage is %s outputDir fileName initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	fileName = sys.argv[2]
	initFrame = int(sys.argv[3])

	msd = PolyMSD(outputDir, fileName, initFrame=initFrame)

	if len(sys.argv) == 3:
		msd.Compute()
		msd.Print()
		
	elif len(sys.argv) == 4:
		idxTad = int(sys.argv[3])
	
		msd.ComputeTad(idxTad)
		msd.PrintTad(idxTad)
