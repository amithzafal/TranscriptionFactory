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


class PolyEndtoEnd():

	def __init__(self, outputDir, initFrame):
		self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)
		
		self.endtoendFile = os.path.join(self.reader.outputDir, "polyEndtoEnd.res")

		if os.path.exists(self.endtoendFile):
			print("Files '%s' and '%s' already exist - aborting" % (self.msdHetFile, self.msdHomFile))
			sys.exit()


	def Compute(self):
		self.endtoend=[]
		
		for i in range(self.reader.N):
			data = next(self.reader)
					
			diff=data.polyPos[0]-data.polyPos[-1]
			self.endtoend.append(np.dot(diff,diff.T)**0.5)
	


	
	def Print(self):
		np.savetxt(self.endtoendFile, self.endtoend)
			


	
if __name__ == "__main__":
	if len(sys.argv) not in [3]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])

	endtoend = PolyEndtoEnd(outputDir, initFrame=initFrame)

	if len(sys.argv) == 3:
		endtoend.Compute()
		endtoend.Print()
		

