#!/bin/bash
#$ -N tissueAnno

cd /hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer

/hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/tools/anaconda/envs/r_env/bin/Rscript ./src/annotation.R
