#!/bin/bash

#$ -N mixer
#$ -o /exports/eddie/scratch/s2093301/coloc/resources/MIXER/out/
#$ -e /exports/eddie/scratch/s2093301/coloc/resources/MIXER/err/
#$ -l h_vmem=16G -pe sharedmem 16
#$ -l h_rt=24:00:00 
#$ -t 1-20

. /etc/profile.d/modules.sh
module purge
module load singularity

DIR=/exports/eddie/scratch/s2093301/coloc/MIXER
MIXER_REFERENCE_FOLDER=$DIR/mixer
MIXER_SIF=$DIR/mixer.sif
cd $DIR/../resources/MIXER/

SCR_TMP=/exports/eddie/scratch/$USER/tmp
mkdir -p "$SCR_TMP"/{apptainer_cache,apptainer_tmp}
export APPTAINER_CACHEDIR="$SCR_TMP/apptainer_cache"
export APPTAINER_TMPDIR="$SCR_TMP/apptainer_tmp"
export APPTAINERENV_TMPDIR="$SCR_TMP"
export APPTAINER_ENABLE_OVERLAY=no
export APPTAINER_NO_MOUNT=tmp
export SINGULARITY_CACHEDIR="$SCR_TMP/apptainer_cache"
export SINGULARITY_TMPDIR="$SCR_TMP/apptainer_tmp"
export SINGULARITY_WORKDIR="$SCR_TMP"
export SINGULARITY_NO_MOUNT=tmp
export TMPDIR="$SCR_TMP"

export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
export OPENBLAS_NUM_THREADS=16

export SINGULARITY_BIND=${MIXER_REFERENCE_FOLDER}/reference:/REF
export MIXER_SIF=$DIR/mixer.sif

export MIXER_COMMON_ARGS="--ld-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld --bim-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim --threads 16"
export REP="${SGE_TASK_ID}"
export EXTRACT="--extract /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep$REP.snps"

export PYTHON="singularity exec --home $PWD:/home $MIXER_SIF python"
export MIXER_PY="$PYTHON /tools/mixer/precimed/mixer.py"

$MIXER_PY test2 $MIXER_COMMON_ARGS --trait1-file MDD_MiXeR.sumstats.gz --trait2-file MCP_MiXeR.sumstats.gz --load-params MDD_vs_MCP.fit.$REP.json --out MDD_vs_MCP.test.$REP

