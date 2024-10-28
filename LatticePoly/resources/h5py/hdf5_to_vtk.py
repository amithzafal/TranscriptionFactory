import os
import sys

import numpy as np

import h5py
from hdf5Reader import hdf5Reader

from vtk import vtkXMLPolyDataWriter
from vtk import vtkPoints
from vtk import vtkFloatArray
from vtk import vtkPolyData
from vtk import vtkCellArray
from vtk import vtkIntArray
from vtk import vtkLine
from vtk import vtkCubeSource

class Hdf5_to_vtk():
        def __init__(self, inputDir, fileName, outputDir, initFrame=-1):
                
                self.reader = hdf5Reader(inputDir, fileName, initFrame, readLiq=True, readPoly=True, backInBox=False)

                os.chdir(outputDir)
                self.outputPath = fileName.split(".")[0] + "_vtk"
                os.makedirs(self.outputPath, exist_ok=True)
        
        def Print(self):
                self.PrintBox()

                for i in range(self.reader.N):
                        data = next(self.reader)
                        self.PrintLiqFrame(data, i)
                        self.PrintPolyFrame(data, i)
                                                                        
                        if (i+1) % 100 == 0:
                                print("Printed %d out of %d configurations" % (i+1, self.reader.N))


        def PrintLiqFrame(self, data, i):

                fileLiqName = 'liq{:05d}.vtp'.format(i+self.reader.initFrame)
                fileLiqPath = os.path.join(self.outputPath,fileLiqName)
                
                points = vtkPoints()
                liqDensity = vtkFloatArray()
                liqDisplacement = vtkFloatArray()
                
                liqDensity.SetName("Density")
                liqDensity.SetNumberOfComponents(1)
                
                liqDisplacement.SetName("Displacement")
                liqDisplacement.SetNumberOfComponents(3)
                        
                for j in range(self.reader.nLiq):
                
                        aveDensity = data.liqDens[j]

                        x = data.liqPos[j][0]
                        y = data.liqPos[j][1]
                        z = data.liqPos[j][2]
                        
                        dx = data.liqDisp[j][0]
                        dy = data.liqDisp[j][1]
                        dz = data.liqDisp[j][2]
                                        
                        points.InsertNextPoint(x, y, z)
                
                        liqDensity.InsertNextValue(aveDensity)
                        liqDisplacement.InsertNextTuple3(dx, dy, dz)
                
                
                polyData = vtkPolyData()
                writer = vtkXMLPolyDataWriter()

                polyData.SetPoints(points)
                
                polyData.GetPointData().AddArray(liqDensity)
                polyData.GetPointData().AddArray(liqDisplacement)

                writer.SetFileName(fileLiqPath)
                writer.SetInputData(polyData)
                
                writer.Write()

        def PrintPolyFrame(self, data, i):

                filePolyName = 'poly{:05d}.vtp'.format(i+self.reader.initFrame)
                filePolyPath = os.path.join(self.outputPath,filePolyName)

                
                points = vtkPoints()
                lines = vtkCellArray()
                
                types = vtkFloatArray()
                painters = vtkFloatArray()
                sisterIDs = vtkIntArray()

                types.SetName("TAD type")
                types.SetNumberOfComponents(1)

                painters.SetName("Painter status")
                painters.SetNumberOfComponents(1)
                
                sisterIDs.SetName("Sister ID")
                sisterIDs.SetNumberOfComponents(1)
                

                for t in range(self.reader.nTad):
                
                        type = data.polyType[t]
                        id = t
                        
                        painter = data.polyPainter[t]
                        
                        points.InsertNextPoint(data.polyPos[t][0], data.polyPos[t][1], data.polyPos[t][2])
                        
                        types.InsertNextValue(type)
                        painters.InsertNextValue(painter)
                        sisterIDs.InsertNextValue(id)
                
                
                for j in range(self.reader.nTad-1):
                        line = vtkLine()
                        
                        line.GetPointIds().SetId(0, j)
                        line.GetPointIds().SetId(1, j+1)
                
                        lines.InsertNextCell(line)

                
                polyData = vtkPolyData()
                writer = vtkXMLPolyDataWriter()

                polyData.SetPoints(points)
                polyData.SetLines(lines)
                
                polyData.GetPointData().AddArray(types)
                polyData.GetPointData().AddArray(painters)
                polyData.GetPointData().AddArray(sisterIDs)

                writer.SetFileName(filePolyPath)
                writer.SetInputData(polyData)
                
                writer.Write()
	        
        def PrintBox(self):
                fileBoxName = 'box.vtp'
                fileBoxPath = os.path.join(self.outputPath,fileBoxName)

                cubeSource = vtkCubeSource()
                
                L = self.reader.boxDim[0]

                cubeSource.SetCenter((L-0.5)/2., (L-0.5)/2., (L-0.5)/2.)
                
                cubeSource.SetXLength(L+0.5)
                cubeSource.SetYLength(L+0.5)
                cubeSource.SetZLength(L+0.5)
                
                cubeSource.Update()

                writer = vtkXMLPolyDataWriter()
                
                writer.SetFileName(fileBoxPath)
                writer.SetInputConnection(cubeSource.GetOutputPort())
                
                writer.Write()

if __name__ == "__main__":
        if len(sys.argv) != 5:
                print("\033[1;31mUsage is %s inputDir fileName outputDir initFrame\033[0m" % sys.argv[0])
                sys.exit()

        inputDir = sys.argv[1]
        fileName = sys.argv[2]
        outputDir = sys.argv[3]
        initFrame = int(sys.argv[4])


        vtk_writer = Hdf5_to_vtk(inputDir, fileName, outputDir, initFrame=initFrame)

        vtk_writer.Print()