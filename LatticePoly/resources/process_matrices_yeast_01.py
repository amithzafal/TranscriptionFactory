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
import math
import numpy as np
from iced import normalization

from utils import msdFFT
from vtkReader import vtkReader
import networkx as nx
import time
import json
from numpy import nanmean
import numba
from numba import jit, int32
from cooler.create import ArrayLoader
import cooler
import warnings
warnings.filterwarnings('ignore')
import pandas as pd


outputDir = sys.argv[1]
name = sys.argv[2]




@numba.jit()
def merge_matrices(outputDir,name):
	matrices=[]
	time=0
	for folder in os.listdir(outputDir)[:]:
		if(folder.endswith('cool')==False and folder.endswith('.cool')==False and folder.endswith('.gz')==False and folder.endswith('.res')==False and folder.startswith(".n")==False):
			print(folder)
			for file_name in os.listdir(outputDir+'/'+folder):
				check=0
				if file_name==name:
					#1print("file name = "+file_name)
					file_path = os.path.join(outputDir+'/'+folder, file_name)
					clr=cooler.Cooler(file_path)
					matr= (np.array(clr.matrix(balance=False)[:]))
					if(len(matr)>0):
						if(time==0):
							matr_sum=matr
						else:
							matr_sum=matr_sum+matr
						#matrices.append(matr)
						time+=1
						check=1
					break;
			if(check==0):
				print("missing "+ name+" from "+folder)





	#rawdata=np.nansum(matrices,axis=0)

	return [matr_sum,clr]
	
	
mtr=merge_matrices(outputDir,name)

chromsizes= mtr[1].chromsizes
bins = cooler.binnify(chromsizes, 1000)
pixels = ArrayLoader(bins, mtr[0], chunksize=10000000)
cooler.create_cooler(outputDir+'/01_'+name,bins,pixels)





