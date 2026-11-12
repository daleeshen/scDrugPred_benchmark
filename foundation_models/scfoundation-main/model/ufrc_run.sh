#!/bin/bash
#SBATCH --job-name=uce-ebd
#SBATCH --output=benchmarking_%j.out
#SBATCH --error=benchmarking_%j.err
#SBATCH --mail-user=wang.qing@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=gpu
#SBATCH --gpus=a100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-socket=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16gb
#SBATCH --distribution=block:block
#SBATCH --time=5:00:00

pwd; hostname; date

miniconda3=/home/wang.qing/miniconda3/bin
export PATH=$miniconda3:$PATH

#export PATH="/home/wang.qing/miniconda3/bin:$PATH

#working directory
#cd /blue/qsong1/wang.qing/benchmark_scLLM/scFoundation-main/scFoundation-main/model/

python dataset_making.py
