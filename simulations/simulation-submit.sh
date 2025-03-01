#!/bin/sh

HOME="/home/students/mnorwood"

SGEOUT="${HOME}/biost571/sgeout/"
OUTPUT="${HOME}/biost571/project3/"
SHELL="${HOME}/repos/biost571-final/simulations/shell.sh"
SCRIPT="${HOME}/repos/biost571-final/simulations/simulation-worker.R"

CONSTANTS="-cwd -N sim571 -j y -o ${SGEOUT} -pe smp 1 -q normal.q ${SHELL} ${SCRIPT}"

# This will do the number of simulations
# passed to this script as the first arg
for i in $(seq 1 1 $1)
do
  ARGS="${i} ${OUTPUT}"
  qsub ${CONSTANTS} ${ARGS}
done
