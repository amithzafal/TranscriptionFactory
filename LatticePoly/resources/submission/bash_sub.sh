#!/bin/bash


##
##  slurm.sh
##  LatticePoly
## 
##  Created by mtortora on 20/06/2022
##  Copyright © 2022 ENS Lyon. All rights reserved.
##

#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

# Relative path to code root directory
ROOTDIR=${SCRIPTDIR}/../..

#!/bin/bash


##
##  slurm.sh
##  LatticePoly
## 
##  Created by mtortora on 20/06/2022
##  Copyright © 2022 ENS Lyon. All rights reserved.
##

#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

# Relative path to code root directory
ROOTDIR=${SCRIPTDIR}/../..

#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/Drosophila/LatticePoly/LatticePoly/resources/MonomerDist_HiC_S_cool.py /scratch/Bio/ddasaro/LatticeData/null_mesc_3D/  1500 3 100
# Set working directory to root
#rm -rf /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/*/*/r_1*
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_120_d0.05 5010 1
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.05 5010 1
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.075 5010 1
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.2 5010 1
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.3 5010 1
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.4 5010 1
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.5 5010 1
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py  /Xnfs/lbmcdb/Jost_team/ddasaro/year3/September23/mitotic_chr4_2/N_70_d0.6 5010 1
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_all_frames_general.py /scratch/Bio/ddasaro/LatticeData/Tads_2L_0.0 500 2
 
#python LatticePoly/LatticePoly/merging_scripts/merge_single_script.py /Xnfs/physbiochrom/ddasaro/year2/January23/segregation_25_long/double/ polyHomMSD

#python Single_cell_scratch.py
#bash /home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/clusters.py /Xnfs/physbiochrom/ddasaro/year3/October23/cohesion_S/null_model_Jsister_0 500
#cp -r  /home/ddasaro/LatticeData/mrt1 /Xnfs/physbiochrom/ddasaro/year3/March24/
#rm -rf /home/ddasaro/LatticeData/S_phase*


# python /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_yeast.py /home/ddasaro/LatticeData/G1_800nm/ r_2.0_0_100full_genome.cool
# python /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_yeast.py /home/ddasaro/LatticeData/G1_800nm/ r_2.0_100_100full_genome.cool
 
#python /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_yeast.py /home/ddasaro/LatticeData/G1_800nm/ r_2.0_200_100full_genome.cool
# python /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_yeast.py /home/ddasaro/LatticeData/G1_800nm/ r_2.0_300_100full_genome.cool
 #python /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_yeast.py /home/ddasaro/LatticeData/G1_800nm/ r_2.0_500_100full_genome.cool

# python /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_yeast.py /home/ddasaro/LatticeData/G1_800nm/ r_2.0_700_100full_genome.cool
# python /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_yeast.py /home/ddasaro/LatticeData/G1_800nm/ r_2.0_900_100full_genome.cool

 #python /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_yeast.py /home/ddasaro/LatticeData/G1_800nm/ r_2.0_500_100full_genome.cool
 #python /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_yeast.py /home/ddasaro/LatticeData/G1_800nm/ r_4.0_500_100full_genome.cool
 #python /home/ddasaro/LatticePoly/LatticePoly/resources/process_matrices_yeast.py /home/ddasaro/LatticeData/G1_800nm/ r_1.0_500_100full_genome.cool

/home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py /Xnfs/abc/nf_scratch/ddasaro/LatticeData/N_150_d_0.06/ 550 5
/home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py /Xnfs/abc/nf_scratch/ddasaro/LatticeData/N_150_d_0.06/ 550 6

/home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py /Xnfs/abc/nf_scratch/ddasaro/LatticeData/N_70_d_0.06/ 550 5
/home/ddasaro/LatticePoly/LatticePoly/resources/process.sh /home/ddasaro/LatticePoly/LatticePoly/resources/MonomerDist_HiC_G1_chr4.py /Xnfs/abc/nf_scratch/ddasaro/LatticeData/N_70_d_0.06/ 550 6
