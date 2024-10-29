#!/bin/bash
##
##  slurm_sweep_tmp.sh
##  LatticePoly
##
##  Created by ppuel on 29/10/2024
##  Copyright © 2022 ENS Lyon. All rights reserved.
##

#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

# Output directory
OUTPUTDIR=Sim1

# Temporary directory
TEMPORARYDIR=/home/${LOGNAME}/Documents/Simulation/tmp/tmp

# Associated scratch directory
SCRATCHDIR=/home/${LOGNAME}/Documents/Simulation/tmp/data

# Script (relative) path
SCRIPTDIR=/home/paulswann/Documents/Simulation/LatticePoly/LatticePoly/resources/submission

# Data directory
XNFSDIR=/home/${LOGNAME}/Documents/Simulation/tmp/archive

# Relative path to code root directory
ROOTDIR=${SCRIPTDIR}/../..

# Set working directory to root
cd ${ROOTDIR}

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTDIR/lib
PYTHONPATH=$PYTHONPATH:$ROOTDIR/resources/h5py
export LD_LIBRARY_PATH
export PYTHONPATH

# Executable path
EXEC=bin/lat

#Values of the parameters
SLURM_ARRAY_TASK_ID=47

JLL_ARRAY=(0.5 0.6 0.7 0.8 0.5 0.6 0.7 0.8 0.5 0.6 0.7 0.8 0.5 0.6 0.7 0.8 0.5 0.6 0.7 0.8 0.5 0.6 0.7 0.8 0.5 0.6 0.7 0.8 0.5 0.6 0.7 0.8 0.5 0.6 0.7 0.8 0.5 0.6 0.7 0.8 0.5 0.6 0.7 0.8 0.5 0.6 0.7 0.8)
JLL=${JLL_ARRAY[$[$SLURM_ARRAY_TASK_ID-1]]}

JLP=0.0
JLPP=0.0
EV=0.0

JLL_VALENCY_ARRAY=(1 1 1 1 2 2 2 2 3 3 3 3 1 1 1 1 2 2 2 2 3 3 3 3 1 1 1 1 2 2 2 2 3 3 3 3 1 1 1 1 2 2 2 2 3 3 3 3)
JLL_VALENCY=${JLL_VALENCY_ARRAY[$[$SLURM_ARRAY_TASK_ID-1]]}

JLP_VALENCY=0
JPL_VALENCY=0
JPPL_VALENCY=0
JLPP_VALENCY=0

LDENS_ARRAY=(0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.02 0.03 0.03 0.03 0.03 0.03 0.03 0.03 0.03 0.03 0.03 0.03 0.03 0.04 0.04 0.04 0.04 0.04 0.04 0.04 0.04 0.04 0.04 0.04 0.04)
LDENS=${LDENS_ARRAY[$[$SLURM_ARRAY_TASK_ID-1]]}

# Temporary directory on tmp
TMPDIR=${TEMPORARYDIR}/${OUTPUTDIR}_JLL_${JLL}_JLL_VALENCY_${JLL_VALENCY}_LDENS_${LDENS}_N_0


# H5 File path
FILEPATH=${TMPDIR}/traj.h5
# Create temporary directory if necessary
[ ! -d "${TMPDIR}" ] && mkdir -p ${TMPDIR}

# Scratch directory on scratch
SCRDIR=${SCRATCHDIR}/${OUTPUTDIR}/JLL/${JLL}/JLL_VALENCY/${JLL_VALENCY}/LDENS/${LDENS}/N/0

# Create scratch directory if necessary
[ ! -d "${SCRDIR}" ] && mkdir -p ${SCRDIR}

# Data directory on Xnfs
DATDIR=${XNFSDIR}/${OUTPUTDIR}/JLL/${JLL}/JLL_VALENCY/${JLL_VALENCY}/LDENS/${LDENS}/N/0

# Create data directory if necessary
[ ! -d "${DATDIR}" ] && mkdir -p ${DATDIR}

# Substitution strings
DIRSUB="s|\(outputDir[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${TMPDIR} ;|;"
FILSUB="s|\(H5filePath[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${FILEPATH} ;|;"
JLLSUB="s|\(Jll[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${JLL} ;|;"
JLPSUB="s|\(Jlp[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${JLP} ;|;"
JLPPSUB="s|\(Jlpp[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${JLPP} ;|;"
EVSUB="s|\(EV[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${EV} ;|;"
JLL_VALENCYSUB="s|\(Jll_Valency[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${JLL_VALENCY} ;|;"
JLP_VALENCYSUB="s|\(Jlp_Valency[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${JLP_VALENCY} ;|;"
JPL_VALENCYSUB="s|\(Jpl_Valency[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${JPL_VALENCY} ;|;"
JPPL_VALENCYSUB="s|\(Jppl_Valency[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${JPPL_VALENCY} ;|;"
JLPP_VALENCYSUB="s|\(Jlpp_Valency[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${JLPP_VALENCY} ;|;"
LDENSSUB="s|\(Ldens[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${LDENS} ;|;"

# Copy input configuration file to output directory, substituting paths and parameter values
sed -e "${DIRSUB}""${FILSUB}""${JLLSUB}""${JLPSUB}""${JLPPSUB}""${EVSUB}""${JLL_VALENCYSUB}""${JLP_VALENCYSUB}""${JPL_VALENCYSUB}""${JPPL_VALENCYSUB}""${JLPP_VALENCYSUB}""${LDENSSUB}" < data/input.cfg > ${TMPDIR}/input.cfg

# Run
./${EXEC} ${TMPDIR}/input.cfg > ${TMPDIR}/log.out

# Perform post-processing analyses
.venv/bin/python3 resources/h5py/LiqDensity.py ${TMPDIR} traj.h5 -1 >> ${TMPDIR}/process.out
.venv/bin/python3 resources/h5py/LiqCluster.py ${TMPDIR} traj.h5 -1 >> ${TMPDIR}/process.out
.venv/bin/python3 resources/h5py/LiqMSD.py ${TMPDIR} traj.h5 -1 >> ${TMPDIR}/process.out
.venv/bin/python3 resources/h5py/LiqPolyCoincidence.py ${TMPDIR} traj.h5 -1 >> ${TMPDIR}/process.out
.venv/bin/python3 resources/h5py/PolyGyration.py ${TMPDIR} traj.h5 -1 >> ${TMPDIR}/process.out
.venv/bin/python3 resources/h5py/PolyMSD.py ${TMPDIR} traj.h5 -1 >> ${TMPDIR}/process.out


# Move processed output files to scratch directory
cp ${TMPDIR}/input.cfg ${SCRDIR}
cp ${TMPDIR}/log.out ${SCRDIR}
cp ${TMPDIR}/process.out ${SCRDIR}
cp ${TMPDIR}/process.h5 ${SCRDIR}

# Move all files to XNFS directory
mv ${FILEPATH} ${DATDIR}/
mv ${TMPDIR}/input.cfg ${DATDIR}
mv ${TMPDIR}/log.out ${DATDIR}
mv ${TMPDIR}/process.out ${DATDIR}
mv ${TMPDIR}/process.h5 ${DATDIR}

# Clean scratch
rm -rf ${TMPDIR}
