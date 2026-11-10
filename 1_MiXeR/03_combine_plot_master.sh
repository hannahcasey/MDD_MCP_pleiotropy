#!/bin/bash
#$ -N bivariate_mixer_figures
#$ -l h_rss=8G
#$ -l h_rt=1:00:00
#$ -l h_rss=8G
#$ -m baes
#$ -M sh.casey@sms.ed.ac.uk
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
export MIXER_FIG="$PYTHON /tools/mixer/precimed/mixer_figures.py"


$MIXER_FIG combine --json ${RESULTS_DIR}/${TRAIT1}.fit.@.json  --out ${RESULTS_DIR}/${TRAIT1}.fit
$MIXER_FIG combine --json ${RESULTS_DIR}/${TRAIT1}.test.@.json  --out ${RESULTS_DIR}/${TRAIT1}.test
$MIXER_FIG combine --json ${RESULTS_DIR}/${TRAIT2}.fit.@.json  --out ${RESULTS_DIR}/${TRAIT2}.fit
$MIXER_FIG combine --json ${RESULTS_DIR}/${TRAIT2}.test.@.json  --out ${RESULTS_DIR}/${TRAIT2}.test
$MIXER_FIG combine --json ${RESULTS_DIR}/${TRAIT1}_vs_${TRAIT2}.fit.@.json  --out ${RESULTS_DIR}/${TRAIT1}_vs_${TRAIT2}.fit
$MIXER_FIG combine --json ${RESULTS_DIR}/${TRAIT1}_vs_${TRAIT2}.test.@.json  --out ${RESULTS_DIR}/${TRAIT1}_vs_${TRAIT2}.test


$MIXER_FIG one --json ${RESULTS_DIR}/${TRAIT1}.fit.json ${RESULTS_DIR}/${TRAIT2}.fit.json --out ${RESULTS_DIR}/${TRAIT1}_and_${TRAIT2}.fit --trait1 ${TRAIT1} ${TRAIT2} --statistic mean std --ext png
$MIXER_FIG one --json ${RESULTS_DIR}/${TRAIT1}.test.json ${RESULTS_DIR}/${TRAIT2}.test.json --out ${RESULTS_DIR}/${TRAIT1}_and_${TRAIT2}.test --trait1 ${TRAIT1} ${TRAIT2} --statistic mean std --ext png
$MIXER_FIG two --json-fit ${RESULTS_DIR}/${TRAIT1}_vs_${TRAIT2}.fit.json --json-test ${RESULTS_DIR}/${TRAIT1}_vs_${TRAIT2}.test.json --out ${RESULTS_DIR}/${TRAIT1}_vs_${TRAIT2} --trait1 ${TRAIT1} --trait1-color 1 --trait2 ${TRAIT2} --trait2-color 2 --statistic mean std --ext png
