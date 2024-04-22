# -*- coding: utf-8 -*-
##
# MonomerDist.py
# LatticePoly
##
# Based on mtortora's  code.
# Copyright © 2021 ENS Lyon. All rights reserved.
##

import os
import sys
import pandas as pd
import numpy as np
import cooler
from scipy.spatial import cKDTree
from cooler.create import ArrayLoader
from vtkReader_multi import vtkReader

from scipy.spatial.distance import pdist, squareform


class theta_arm():
	def __init__(self, outputDir, initFrame):
		self.readers=[]
		for i in range(0,17):
			self.readers.append(vtkReader(outputDir,i,initFrame,readLiq=False, readPoly=True))
		self.costhetaFile = os.path.join(outputDir,str(initFrame)+"_"+str(interval)+"_cos_thetha_arm_full_genome.res")
		self.costheta_err_File = os.path.join(outputDir,str(initFrame)+"_"+str(interval)+"_err_cos_thetha_arm_full_genome.res")
		self.centromeres=np.loadtxt("/home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/data/centromeres.in",dtype=int)
		self.costheta_tableFile = os.path.join(outputDir,str(initFrame)+"_"+str(interval)+"_cos_thetha_arm_table_full_genome.res")
		#self.Nchain=0
		#for t in range(self.reader.nTad):
		#	if(self.reader.status[t]==-1 or self.reader.status[t]==0):
		#		self.Nchain+=1
		#self.timepoint=initFrame #find frame where replication reach the desired percentage
		#if(self.reader.nTad>self.Nchain):
		#	print("Chromosome already replicating ")
		#self.timepoint=init_time #find frame where replication starts
		
				

		#print(self.timepoint)
		#self.reader = vtkReader(outputDir, self.timepoint,readLiq=False, readPoly=True)
		#restarted vtk reader from middle frame of desired percentage
		#compute the hic for the minutes
		self.Compute(interval)

		self.Print()
			






		
	#NB here is not the finalFrame but the number of iterations
	def Compute(self,finalFrame):
		self.costheta_table=np.zeros((17,finalFrame))
		for i in range(0, finalFrame):
			self.ProcessFrame(i)
			if (i+1) % 10 == 0:
				print("Processed %d out of %d configurations" %
					  (i+1, finalFrame))


#self.contactProb=np.rint(self.contactProb/finalFrame)

	def ProcessFrame(self, i):
		for chrom,reader in enumerate(self.readers):
			if(chrom!=12):
				#print(self.centromeres)
				data = next(reader)
				SPB=data.polyPos[self.centromeres[chrom]]
				r_cm1=np.nanmean(data.polyPos[:self.centromeres[chrom]],axis=0)-SPB
				r_cm2=np.nanmean(data.polyPos[self.centromeres[chrom]:],axis=0)-SPB
				self.costheta_table[chrom,i]=np.arccos(np.dot(r_cm1,r_cm2)/((np.sqrt(np.dot(r_cm1,r_cm1))*np.sqrt(np.dot(r_cm2,r_cm2)))))


			

					



	def Print(self):
		np.savetxt(self.costheta_tableFile,self.costheta_table)
		np.savetxt(self.costhetaFile,np.array(np.nanmean(self.costheta_table,axis=1)).T)
		np.savetxt(self.costheta_err_File,np.array(np.std(self.costheta_table,axis=1)/interval**0.5).T)



if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("\033[1;31mUsage is %s outputDir initFrame interval \033[0m" % sys.argv[0])
		sys.exit()
	
	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	interval =	int(sys.argv[3])
	#init_time=int(sys.argv[4])


theta = theta_arm(outputDir, initFrame=initFrame)

