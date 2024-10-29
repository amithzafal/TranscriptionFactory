##
##  LiqDensity.py
##  LatticePoly
##
##  Created by ppuel on 18/10/2024.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np

from hdf5Reader import hdf5Reader
import h5py

class LiqDensity():

	def __init__(self, outputDir, fileName, initFrame, threshold=0.5):
		
		self.reader = hdf5Reader(outputDir, fileName, initFrame=-1, readLiq=True, readPoly=False, backInBox=False)
		self.filePath = os.path.join(outputDir, "process.h5")
		
		self.threshold = threshold


	def Compute(self):
		self.meanHist = np.zeros(self.reader.N, dtype=np.float32)
		self.stdHist = np.zeros(self.reader.N, dtype=np.float32)

		for i in range(self.reader.N):
			self.ProcessFrame(i)
									
			if (i+1) % 100 == 0:
				print("Processed %d out of %d configurations" % (i+1, self.reader.N))

			
	def ProcessFrame(self, i):
		data = next(self.reader)
		meastdDensnDens = data.liqDens.sum()
		stdDens = np.square(data.liqDens - data.liqDens.mean()).sum()
		self.meanHist[i] = np.count_nonzero(data.liqDens > self.threshold)
		self.stdHist[i] = stdDens

	
	def Print(self):

		self.reader.Close()
		file = h5py.File(self.filePath,'a')
		
		file.create_dataset("liqMean", data = self.meanHist / self.reader.nLiq)
		file.create_dataset("liqSTD", data = np.sqrt(self.stdHist / self.reader.nLiq))
		print("\033[1;32mPrinted liquid mean densities to '%s'\033[0m" % "liqMean")
		print("\033[1;32mPrinted liquid density STDs fraction to '%s'\033[0m" % "liqSTD")

		file.close()

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("\033[1;31mUsage is %s outputDir fileName initFrame\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	fileName = sys.argv[2]
	initFrame = int(sys.argv[3])

	density = LiqDensity(outputDir, fileName, initFrame=initFrame)

	density.Compute()
	density.Print()
