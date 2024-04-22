#!/bin/bash

##
##  slurm.sh
##  LatticePoly
##
##  Created by ddasaro on the model of mtortora on 20/06/2022
##  Copyright © 2022 ENS Lyon. All rights reserved.
##


#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

# Relative path to code root directory
ROOTDIR=${SCRIPTDIR}/../..
echo "The present scriptdir directory is ${SCRIPTDIR}";



# Set working directory to root
cd ${ROOTDIR}

# Executable path
EXEC=bin/lat

# Data directory on local disk
DATDIR=data/output

# Name of output directory (based on job name & task ID)
ID=$(echo ${SLURM_ARRAY_TASK_ID} | awk '{printf("%04d\n", $1-1)}')

OUTDIR=${SLURM_JOB_NAME}${ID}

# Temporary directory on scratch
[ ! -d "${SCRATCHDIR}/${LOGNAME}/LatticeData/${SLURM_JOB_NAME}/" ] && mkdir -p ${SCRATCHDIR}/${LOGNAME}/LatticeData/${SLURM_JOB_NAME}/
mkdir ${SCRATCHDIR}/${LOGNAME}/LatticeData/${SLURM_JOB_NAME}/${OUTDIR}
TMPDIR=${SCRATCHDIR}/${LOGNAME}/LatticeData/${SLURM_JOB_NAME}/${OUTDIR}

# Create output directory if necessary
#[ ! -d "${TMPDIR}" ] && mkdir -p ${TMPDIR}

# Substitution strings
DIRSUB="s|\(outputDir[[:space:]]*=[[:space:]]*\)\(.*;\)|\1${TMPDIR} ;|;"

# Copy input configuration file to output directory, substituting paths and parameter values
sed -e "${DIRSUB}" < /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/submission/${TEMP}/input.cfg > ${TMPDIR}/input.cfg


# Run
./${EXEC} ${TMPDIR}/input.cfg > ${TMPDIR}/log.out

#Analysis


#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  0 2 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  100 2 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  200 2100

#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  300 2 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  400 2 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  500 1 400
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  500 2 400
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  500 3 400
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  500 4 400
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  500 2 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  700 2 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome.py  ${TMPDIR}  900 2 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  500 2 100
python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  500 3 100
python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  600 3 100
python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  700 3 100
python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  900 3 100
python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  1100 3 100

#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  500 4 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  500 5 100

#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  500 2 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  500 3 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  500 4 100
#python3 /home/ddasaro/Yeast_full_genome/LatticePoly/LatticePoly/resources/MonomerDist_HiC_full_genome_S_phase.py ${TMPDIR}  500 5 100




# Move SLURM output files to data directory
mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out ${TMPDIR}
mv ${SLURM_SUBMIT_DIR}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err ${TMPDIR}

# Create data directory on local disk
#[ ! -d "${DATDIR}" ] && mkdir -p ${DATDIR}

# Archive output files to home directory
tar --transform "s|^|${OUTDIR}/|" -czf /Xnfs/physbiochrom/ddasaro/year3/April24/S_phase_full_genome/${OUTDIR}.tar.gz -C ${TMPDIR} .
#tar -xzf ${DATDIR}/${OUTDIR}.tar.gz -C /Xnfs/lbmcdb/Jost_team/ddasaro/year2/April23/Ring/100
#rm -rf {DATDIR}/${OUTDIR}.tar.gz

# Clean scratch
#rm -rf ${TMPDIR}


