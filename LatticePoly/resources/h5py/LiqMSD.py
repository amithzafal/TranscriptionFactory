##
##  LiqMSD.py
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


class LiqMSD():

	def __init__(self, outputDir, fileName, initFrame):
		self.reader = hdf5Reader(outputDir, fileName, initFrame=-1, readLiq=True, readPoly=False, backInBox=False)
		self.filePath = os.path.join(outputDir, "process.h5")


	def Compute(self):
		vMem = psutil.virtual_memory()
		sizeTot = self.reader.N * self.reader.liqPos.nbytes

		if sizeTot < vMem.available:
			self.cumulDist = 0
			self.liqPosInit = self.reader.liqPos

			posHist = self.ReadHist()

			for idxSpin in range(self.reader.nLiq):
				self.cumulDist += msdFFT(posHist[:, idxSpin])

				if (idxSpin+1) % 1000 == 0:
					print("Processed %d out of %d spins" % (idxSpin+1, self.reader.nLiq))

		else:
			print("Memory overflow likely - reduce chosen number of frames")
			sys.exit()


	def ReadHist(self):
		posHist = np.zeros((self.reader.N, self.reader.nLiq, 3), dtype=np.float32)

		for i in range(self.reader.N):
			data = next(self.reader)
			posHist[i] = self.liqPosInit + data.liqDisp

		return posHist


	def Print(self):
		msdLiq = self.cumulDist / self.reader.nLiq

		self.reader.Close()
		file = h5py.File(self.filePath, 'a')
		file.create_dataset("liqMSD", data = msdLiq)
	
		print("\033[1;32mPrinted liquid MSDs to '%s'\033[0m" % "liqMSD")
		
		file.close()


if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("\033[1;31mUsage is %s outputDir fileName initFrame\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	fileName = sys.argv[2]
	initFrame = int(sys.argv[3])

	msd = LiqMSD(outputDir, fileName, initFrame=initFrame)

	msd.Compute()
	msd.Print()
