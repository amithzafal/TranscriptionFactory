##
##  debug.py
##  LatticePoly
##
##  Created by ppuel on 23/10/2024.
##  Copyright © 2024 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np
from vtkReader import vtkReader
from hdf5Reader import hdf5Reader
import h5py as h5
import matplotlib.pyplot as plt

cwd = os.getcwd()

os.chdir("data/output/testVTK")

list_file = os.listdir("test2")

Reader_from_vtk = vtkReader("test2",-1,True,True,True)
Reader_from_hdf5 = hdf5Reader("","test.h5",-1,True,True,True)

# for i in range(10):
#         vtk_data = next(Reader_from_vtk)
#         hdf5_data = next(Reader_from_hdf5)
#         for j in range(len(vtk_data.polyPos)):
#                 print(vtk_data.polyPos[j],hdf5_data.polyPos[j],np.sum(np.square(vtk_data.polyPos[j]-hdf5_data.polyPos[j])))

#         print()

for file in list_file:
        if file.split(".")[-1] == "res":
                print("\n********************\n")
                print(file)
                VTKdata = np.loadtxt(os.path.join("test2",file))

                h5file = h5.File("test.h5","r")
                h5data = np.array(h5file[file.split('.')[0]])

                h5file.close()

                square_error = np.sum(np.square((VTKdata-h5data)))
                print(file,square_error)
                if square_error > 1e-10:
                        for i in range(len(h5data)):
                                print(VTKdata[i],h5data[i])



# ax = plt.figure().add_subplot(projection='3d')
#         ax.scatter(vtk_data.polyPos[:,0], vtk_data.polyPos[:,1], vtk_data.polyPos[:,2], label='hdf5')
#         ax.scatter(hdf5_data.polyPos[:,0], hdf5_data.polyPos[:,1], hdf5_data.polyPos[:,2], label='vtk')
#         ax.legend()
#         ax.set_xlim([0,5])
#         ax.set_ylim([0,5])
#         ax.set_zlim([0,5])
#         plt.show()"