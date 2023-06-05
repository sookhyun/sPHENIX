#!/bin/csh
set workdir = /gpfs/mnt/gpfs02/sphenix/user/shlee/matching/macro 
cd $workdir
echo "Executing run_select.sh in " $workdir

mkdir -p LOGS
mkdir -p OUTPUT

source /opt/sphenix/core/bin/setup_local.csh /sphenix/user/shlee/install

root -l -b -q Fun4All_JetAna.C\(0,$1\)

# 1: file index


