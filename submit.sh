#!/bin/bash
#SBATCH --job-name=submit-imputation-pipeline
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load nextflow
nextflow run main.nf -c nextflow.config