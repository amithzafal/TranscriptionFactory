#!/bin/bash

##
##  slurm.sh
##  LatticePoly
##
##  Created by ppuel on 23/10/2024
##  Copyright © 2024 ENS Lyon. All rights reserved.
##

#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

# Relative path to code root directory
ROOTDIR=~/Documents/Simulation/LatticePoly/LatticePoly/

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${ROOTDIR}/lib
export PYTHONPATH=${ROOTDIR}/lib/vtk:$PYTHONPATH
# Set working directory to root
cd ${ROOTDIR}


# Executable path
EXEC=bin/lat

# Run
./${EXEC} data/input.cfg > data/log.out

${ROOTDIR}/.venv/bin/python3 resources/h5py/LiqCluster.py ${ROOTDIR}/data/output/testVTK test.h5 -1
${ROOTDIR}/.venv/bin/python3 resources/h5py/LiqDensity.py ${ROOTDIR}/data/output/testVTK test.h5 -1
${ROOTDIR}/.venv/bin/python3 resources/h5py/LiqMSD.py ${ROOTDIR}/data/output/testVTK test.h5 -1
${ROOTDIR}/.venv/bin/python3 resources/h5py/LiqPolyCoincidence.py ${ROOTDIR}/data/output/testVTK test.h5 -1
${ROOTDIR}/.venv/bin/python3 resources/h5py/PolyGyration.py ${ROOTDIR}/data/output/testVTK test.h5 -1
${ROOTDIR}/.venv/bin/python3 resources/h5py/PolyMSD.py ${ROOTDIR}/data/output/testVTK test.h5 -1

python3 resources/LiqCluster.py ${ROOTDIR}/data/output/testVTK/test2 -1
python3 resources/LiqDensity.py ${ROOTDIR}/data/output/testVTK/test2 -1
python3 resources/LiqMSD.py ${ROOTDIR}/data/output/testVTK/test2 -1
python3 resources/LiqPolyCoincidence.py ${ROOTDIR}/data/output/testVTK/test2 -1
python3 resources/PolyGyration.py ${ROOTDIR}/data/output/testVTK/test2 -1
python3 resources/PolyMSD.py ${ROOTDIR}/data/output/testVTK/test2 -1