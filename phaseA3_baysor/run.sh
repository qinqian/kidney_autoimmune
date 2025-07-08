#!/bin/bash -ex

#SBATCH --job-name=baysor1
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:0:0
#SBATCH --mem=256G
#SBATCH --output=log%J.out
#SBATCH --error=log%J.err

cd /data/wei/qq06/xenium/phaseA3_baysor/../data/kidney/20241025__200743__BWH_20241025_SHRUTI_RACHEL/output-XETG00392__0045655__BS21-N65682A2__20241025__201009/

JULIA_NUM_THREADS=12 baysor run transcripts.parquet -m 20  -x x_location -y y_location -z z_location -g feature_name --n-clusters=4 -s 20 --scale-std=25% 

