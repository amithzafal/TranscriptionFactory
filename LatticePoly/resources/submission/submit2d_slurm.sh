#!/bin/bash

##
##  submit_slurm.sh
##  LatticePoly
##
##  Created by mtortora on 20/06/2022.
##  Copyright © 2022 ENS Lyon. All rights reserved.
##

# Max. walltime
WTIME=8-00:00:00

# Partition
QUEUE="Cascade"

# Max. memory per task
MAXMEM="1G"

# Associated scratch directory
SCRATCHDIR=/home

# Script (relative) path
SCRIPTDIR=$(dirname "$0")

# sbatch arguments
QARGS="--ntasks=1 --mem=${MAXMEM} -a 1-$4 -t ${WTIME} -p ${QUEUE}"

# sbatch variables
QVARS="PARAM=$1,MIN_VAL=$2,MAX_VAL=$3,SCRATCHDIR=${SCRATCHDIR},SCRIPTDIR=${SCRIPTDIR}"

# Check input parameters and submit
if [ "$#" -eq "8" ]; then
	if [ $( expr $4 % 48 ) -eq "0" ]; then
		for i in $(seq $8); do
			VAL2=$(echo $6 $7 1 $8 $i | awk '{printf("%.4f\n", $1+($2-$1)/($4-$3)*($5-$3))}')
			QVARS2="PARAM2=$5,VAL2=${VAL2}"
			
			sbatch ${QARGS} -J $1$5$i --export=${QVARS},${QVARS2} ${SCRIPTDIR}/slurm_sweep.sh
		done
	else
		echo "numJob1 must be multiple of 48 (got $4)"
	fi
else
	echo -e "\033[1;31mUsage is $0 paramName1 minVal1 maxVal1 numJob1 paramName2 minVal2 maxVal2 numJob2\033[0m"
fi
