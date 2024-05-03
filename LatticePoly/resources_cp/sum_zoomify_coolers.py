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
cooler_name = sys.argv[2]




@numba.jit()
def sum_matrices(outputDir,namcooler_name):
	time=0
	processed=0
	for folder in os.listdir(outputDir)[:2]:
		if(folder.endswith('.scool')==False and folder.endswith('.cool')==False and folder.endswith('.gz')==False and folder.endswith('.res')==False and folder.startswith(".n")==False):
			print(folder)
			for file_name in os.listdir(outputDir+'/'+folder):
				check=0
				if file_name==name:
					#1print("file name = "+file_name)
					file_path = os.path.join(outputDir+'/'+folder, file_name)
                    			if(processed==0):
                        			clr_first=cooler.Cooler(file_path)
                        			matrix=clr_first.matrix(balance=False)[:]
						processed+=1
					else:
                        			clr=cooler.Cooler(file_path)
                        			matrix+=clrr.matrix(balance=False)[:]
						break;
				if(check==0):
					print("missing "+ name+" from "+folder)
	return [matrix,clr_first]

mtr=sum_matrices(outputDir,cooler_name)

chromsizes= mtr[1].chromsizes
bins = cooler.binnify(chromsizes, 1000)
pixels = ArrayLoader(bins, mart[0], chunksize=10000000)
cooler.create_cooler(outputDir+'/'+cooler_name,bins,pixels)

cooler.zoomify_cooler(
    outputDir+'/'+cooler_name,
    outputDir+'/'+cooler_name[:-5]+".mcool",
    [1000,2000,4000,8000,16000,32000],
    chunksize=100,
    nproc=10
)




