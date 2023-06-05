#!/bin/csh 

condor_submit condor.job -append "Arguments=$1" -append "Log=/gpfs/mnt/gpfs02/sphenix/user/shlee/matching/macro/LOGS/photonjet_$1.log" -append "Output=/gpfs/mnt/gpfs02/sphenix/user/shlee/matching/macro/LOGS/photonjet_$1.out" -append "Error=/gpfs/mnt/gpfs02/sphenix/user/shlee/matching/macro/LOGS/photonjet_$1.err"
