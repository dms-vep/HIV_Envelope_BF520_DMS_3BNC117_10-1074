#!/bin/bash

#SBATCH -c 1

snakemake --executor slurm --jobs 200 --default-resources -s dms-vep-pipeline-3/Snakefile --use-conda --rerun-incomplete