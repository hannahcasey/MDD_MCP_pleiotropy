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

mkdir -p /exports/eddie/scratch/s2093301/tmp/{apptainer_cache,apptainer_tmp}
export APPTAINER_CACHEDIR=/exports/eddie/scratch/s2093301/tmp/apptainer_cache
export APPTAINER_TMPDIR=/exports/eddie/scratch/s2093301/tmp/apptainer_tmp
export SINGULARITY_CACHEDIR="$APPTAINER_CACHEDIR"
export SINGULARITY_TMPDIR="$APPTAINER_TMPDIR"
export TMPDIR=/exports/eddie/scratch/s2093301/tmp

export SINGULARITY_BIND=${MIXER_REFERENCE_FOLDER}/reference:/REF
export MIXER_SIF=$DIR/mixer.sif

export MIXER_COMMON_ARGS="--ld-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld --bim-file /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim --threads 16"
export REP="${SGE_TASK_ID}"
export EXTRACT="--extract /REF/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep$REP.snps"

export PYTHON="singularity exec --home $PWD:/home $MIXER_SIF python"
export MIXER_PY="$PYTHON /tools/mixer/precimed/mixer.py"


$MIXER_PY fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file MDD_MiXeR.sumstats.gz --out MDD.fit.$REP
$MIXER_PY fit1 $MIXER_COMMON_ARGS $EXTRACT --trait1-file MCP_MiXeR.sumstats.gz --out MCP.fit.$REP
$MIXER_PY fit2 $MIXER_COMMON_ARGS $EXTRACT --trait1-file MDD_MiXeR.sumstats.gz --trait2-file MCP_MiXeR.sumstats.gz --trait1-params MDD.fit.$REP.json --trait2-params MCP.fit.$REP.json --out MDD_vs_MCP.fit.$REP

$MIXER_PY test1 $MIXER_COMMON_ARGS --trait1-file MDD_MiXeR.sumstats.gz --load-params MDD.fit.$REP.json --out MDD.test.$REP
$MIXER_PY test1 $MIXER_COMMON_ARGS --trait1-file MCP_MiXeR.sumstats.gz --load-params MCP.fit.$REP.json --out MCP.test.$REP
$MIXER_PY test2 $MIXER_COMMON_ARGS --trait1-file MDD_MiXeR.sumstats.gz --trait2-file MCP_MiXeR.sumstats.gz --load-params MDD_vs_MCP.fit.$REP.json --out MDD_vs_MCP.test.$REP

