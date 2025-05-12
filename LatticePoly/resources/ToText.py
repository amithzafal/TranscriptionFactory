import os
import sys
import numpy as np
from scipy.spatial import cKDTree
from vtkReader import vtkReader

class ToText:

    def __init__(self, outputDir, chrom1, chrom2, chrom3, initFrame):
        self.reader1 = vtkReader(outputDir, chrom1, initFrame, readLiq=False, readPoly=True)
        self.reader2 = vtkReader(outputDir, chrom2, initFrame, readLiq=False, readPoly=True)
        self.reader3 = vtkReader(outputDir, chrom3, initFrame, readLiq=False, readPoly=True)
	
    def Compute(self):
        for i in range(self.reader1.N):
            self.outputFile = os.path.join(outputDir, f"time_{i}_pos.read_data")
            if os.path.exists(self.outputFile):
                print(f"File '{self.outputFile}' already exists - aborting")
                sys.exit()
            self.ProcessFrame()

    def ProcessFrame(self):
        data1 = next(self.reader1)
        data2 = next(self.reader2)
        data3 = next(self.reader3)
        
        pos1 = data1.polyPos
        pos2 = data2.polyPos
        pos3 = data3.polyPos
        
        coords = np.vstack([pos1, pos2, pos3])
        np.savetxt(self.outputFile, coords)
    def Print(self):
        
        print(f"\033[1;32mPrinted positions to '{self.outputFile}'\033[0m")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("\033[1;31mUsage is %s outputDir chrom1 chrom2 chrom3 initFrame\033[0m" % sys.argv[0])
        sys.exit()

    outputDir = sys.argv[1]
    chrom1 = sys.argv[2]
    chrom2 = sys.argv[3]
    chrom3 = sys.argv[4]
    initFrame = int(sys.argv[5])

    analyzer = ToText(outputDir, chrom1, chrom2, chrom3, initFrame=initFrame)
    analyzer.Compute()
    analyzer.Print()
