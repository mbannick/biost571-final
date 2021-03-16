#!/bin/sh

HOME="/home/students/mnorwood"

SGEOUT="${HOME}/biost571/sgeout/"
OUTPUT="${HOME}/biost571/project/"
SHELL="${HOME}/repos/biost571-final/simulations/shell.sh"
SCRIPT="${HOME}/repos/biost571-final/simulations/simulation-worker.R"

CONSTANTS="-cwd -N sim571 -j y -o ${SGEOUT} -pe smp 1 -q normal.q ${SHELL} ${SCRIPT}"

FITY="exchangeable"
LINK="identity"

ARGS="-corstr ${FITY} -fit_link ${FITY} -directory ${OUTPUT}"

CONSTRAIN="-rho 100"

# This will do the number of simulations
# passed to this script as the first arg
for i in $(seq $1 1 $2)
do
  echo $i
  ARGS="${ARGS} -sim ${i}"
  for N in 100 250
  do
    echo $N
    ARGS="${ARGS} -n_sub ${N}"
    for MAXOBS in 1 5
    do
      echo $MAXOBS
      MINOBS=$((${MAXOBS}+5))
      ARGS="${ARGS} -min_obs ${MINOBS} -max_obs ${MAXOBS}"
      for CORX in independence block
      do
        echo $CORX
        ARGS="${ARGS} -corr_X ${CORX}"
        for ALPHAY in 0.0 0.5
        do
          echo $ALPHAY
          ARGS="${ARGS} -alpha_obs ${ALPHAY}"
          for NPRED in 40 80
          do
            echo $NPRED
            ARGS="${ARGS} -n_pred ${NPRED}"

            qsub ${CONSTANTS} ${ARGS}
            qsub ${CONSTANTS} ${ARGS} ${CONSTRAIN}

            for GAMMA in 1 2
            do
              echo $GAMMA
              PEN="-gamma ${GAMMA}"
              for LAMBDA in 0.5 1 1.5 2
              do
                echo $LAMBDA
                PEN="${PEN} -lambda ${LAMBDA}"

                qsub ${CONSTANTS} ${ARGS} ${PEN}
                qsub ${CONSTANTS} ${ARGS} ${PEN} ${CONSTRAIN}
              done
            done
          done
        done
      done
    done
  done
done
