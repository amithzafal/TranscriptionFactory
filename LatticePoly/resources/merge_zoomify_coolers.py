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
def create_cooler_list(outputDir,namcooler_name):
	matrices=[]
	time=0
	for folder in os.listdir(outputDir)[:]:
		if(folder.endswith('cool')==False and folder.endswith('.cool')==False and folder.endswith('.gz')==False and folder.endswith('.res')==False and folder.startswith(".n")==False):
			print(folder)
			for file_name in os.listdir(outputDir+'/'+folder):
				if file_name==cooler_name:
					#1print("file name = "+file_name)
					file_path = os.path.join(outputDir+'/'+folder, file_name)
					matrices.append(file_path)
					break;
				#print("missing "+ cooler_name +" from "+folder)
	return matrices

coolers_uri=create_cooler_list(outputDir,cooler_name)
print(len(coolers_uri))
cooler.merge_coolers(outputDir+'/'+cooler_name,coolers_uri,mergebuf=100)
print("merging complete")
cooler.zoomify_cooler(outputDir+'/'+cooler_name,outputDir+'/'+cooler_name[:-5]+".mcool",[1000,2000,4000,8000,16000,32000],chunksize=100,nproc=10)




