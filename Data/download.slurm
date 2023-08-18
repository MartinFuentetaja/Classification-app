#!/bin/bash
#SBATCH --job-name=variant_summary_download
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --begin=2023-06-16T09:22:00


./variant_summary_download.sh
