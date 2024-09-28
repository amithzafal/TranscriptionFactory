import os
import sys
import psutil
import math
import numpy as np


outputDir = sys.argv[1]



mrt=[]
for folder in os.listdir(outputDir)[:]:
	current_mrt=[]
	print(folder)
	if(folder.endswith('cool')==False and folder.endswith('.gz')==False and folder.endswith('.res')==False):
		for chrom_n in range(1,18):
			if(np.sum(np.array(os.listdir(outputDir+'/'+folder))=="chrom_"+str(chrom_n-1)+"_timing.res")):
				if 0==0:
					file_path = outputDir+'/'+folder+"/chrom_"+str(chrom_n-1)+"_timing.res"
					if  os.path.getsize(file_path) != 0:
						current_mrt.extend(list(np.loadtxt(file_path)[:]))
		if(len(current_mrt)==12063):						
			mrt.append(current_mrt)
				
print(np.shape(np.array(mrt)))                
print(np.shape(np.sum(np.array(mrt),axis=0)))

np.savetxt(outputDir+"/mrt_raw.res",np.sum(np.array(mrt),axis=0))
