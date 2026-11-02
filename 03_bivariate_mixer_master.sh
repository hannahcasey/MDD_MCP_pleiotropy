#!/bin/bash

#$ -N bivariate_mixer
#$ -l h_rss=16G -pe sharedmem 16
#$ -l h_rt=48:00:00 
#$ -t 1-20
#$ -m baes
#$ -M h.casey@sms.ed.ac.uk
#$ -cwd

. /etc/profile.d/modules.sh

module load apptainer

SCR_TMP=/exports/eddie/scratch/$USER/tmp
mkdir -p "$SCR_TMP"/{apptainer_cache,apptainer_tmp}

export APPTAINER_CACHEDIR="$SCR_TMP/apptainer_cache"
export APPTAINER_TMPDIR="$SCR_TMP/apptainer_tmp"
export MIXER_SIF=/exports/eddie/scratch/s2093301/MDD_MCP/resources/mixer/mixer.sif
export PYTHON="singularity exec --home $PWD:/home $MIXER_SIF python"
export MIXER_PY="$PYTHON /tools/mixer/precimed/mixer.py"

export MIXER_COMMON_ARGS="--ld-file /exports/eddie/scratch/s2093301/MDD_MCP/resources/mixer/1000G/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld --bim-file /exports/eddie/scratch/s2093301/MDD_MCP/resources/mixer/1000G/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim --threads 16"
export REP="${SGE_TASK_ID}"
export EXTRACT="--extract /exports/eddie/scratch/s2093301/MDD_MCP/resources/mixer/1000G/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep$REP.snps"

if [ ${RUN_UNIVARIATE_1} = 1 ]; then
	$MIXER_PY fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file ${SUMSTATS_DIR}/${SUMSTATS_FILE1} --out ${RESULTS_DIR}/${TRAIT1}.fit.$REP
	$MIXER_PY test1 $MIXER_COMMON_ARGS --trait1-file ${SUMSTATS_DIR}/${SUMSTATS_FILE1} --load-params ${RESULTS_DIR}/${TRAIT1}.fit.$REP.json --out ${RESULTS_DIR}/${TRAIT1}.test.$REP
fi

if [ ${RUN_UNIVARIATE_2} = 1 ]; then	
$MIXER_PY fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file ${SUMSTATS_DIR}/${SUMSTATS_FILE2} --out ${RESULTS_DIR}/${TRAIT2}.fit.$REP
$MIXER_PY test1 $MIXER_COMMON_ARGS --trait1-file ${SUMSTATS_DIR}/${SUMSTATS_FILE2} --load-params ${RESULTS_DIR}/${TRAIT2}.fit.$REP.json --out ${RESULTS_DIR}/${TRAIT2}.test.$REP
fi

$MIXER_PY fit2 $MIXER_COMMON_ARGS $EXTRACT --trait1-file ${SUMSTATS_DIR}/${SUMSTATS_FILE1} --trait2-file ${SUMSTATS_DIR}/${SUMSTATS_FILE2} --trait1-params ${RESULTS_DIR}/${TRAIT1}.fit.$REP.json --trait2-params ${RESULTS_DIR}/${TRAIT2}.fit.$REP.json --out ${RESULTS_DIR}/${TRAIT1}_vs_${TRAIT2}.fit.$REP
$MIXER_PY test2 $MIXER_COMMON_ARGS --trait1-file ${SUMSTATS_DIR}/${SUMSTATS_FILE1} --trait2-file ${SUMSTATS_DIR}/${SUMSTATS_FILE2} --load-params ${RESULTS_DIR}/${TRAIT1}_vs_${TRAIT2}.fit.$REP.json --out ${RESULTS_DIR}/${TRAIT1}_vs_${TRAIT2}.test.$REP