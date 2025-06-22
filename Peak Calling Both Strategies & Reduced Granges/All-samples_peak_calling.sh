#!/bin/bash

#SBATCH --account=kellystr_1320
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --time=48:00:00

bash -c "/home1/yulongqi/macs3_env/bin/macs3 callpeak -t /project/kellystr_1320/GEO_Data/*.tsv.gz -g 2.7e+09 -f BED --nomodel --extsize 200 --shift -100 -n SeuratProject --outdir /project/kellystr_1320/GEO_Data/peaks_allsamples"

