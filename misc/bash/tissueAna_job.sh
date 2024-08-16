#!/bin/bash
#$ -N tissueAna
#$ -pe smp 64

cd /hsfscqjf1/ST_CQ/P21Z10200N0096/CRC/lizehua/test/lungcancer

make tissueAna
