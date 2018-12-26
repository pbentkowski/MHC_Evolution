#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=4gb
#SBATCH --time=16:00:00
#SBATCH --partition=standard


module load openmpi/4.0.0_gcc620
module load boost/1.63.0-gcc620_py2712

# Ustawiamy zmienna $TMPDIR
export TMPDIR="/tmp/lustre_shared/${USER}/${SLURM_JOBID}"

# Ustawiamy zmienne aplikacji
export SCR=${TMPDIR}
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Ustawiamy zmienne pomocnicze
#INPUT_DIR="input"
#OUTPUT_DIR="output"
#OUTPUT_FILE="OUTPUT"

# Tworzymy katalog tymczasowy
mkdir -p ${TMPDIR}

# Kopiujemy dane wejsciowe do katalogu wskazywanego zmienna $TMPDIR
cp ${SLURM_SUBMIT_DIR}/* ${TMPDIR}

# Przechodzimy do katalogu $TMPDIR
cd $TMPDIR

# Naglowek

cat << EOF
-------------------------------------------------------------------------------

Start of calculations [$(date)]
EOF

# SEKCJA RUN Wykonujemy obliczenia

