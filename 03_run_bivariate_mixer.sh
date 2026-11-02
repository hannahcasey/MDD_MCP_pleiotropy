export APPTAINER_BIND=/exports/eddie/scratch/s2093301/MDD_MCP/
export MIXER_DIR=/exports/eddie/scratch/s2093301/MDD_MCP/resources/mixer/
export SCRIPT_DIR=/exports/eddie/scratch/s2093301/MDD_MCP/bin
export SUMSTATS_DIR=/exports/eddie/scratch/s2093301/MDD_MCP/resources/sumstats/

# MCP X MDD (EUR meta-analysis - UKB and 23andMe removed)
export SUMSTATS_FILE1=MCP_MiXeR.sumstats.gz
export TRAIT1=MCP
export RUN_UNIVARIATE_1=1

export SUMSTATS_FILE2=MDD_MiXeR.sumstats.gz
export TRAIT2=MDD
export RUN_UNIVARIATE_2=1

export RESULTS_DIR=/exports/eddie/scratch/s2093301/MDD_MCP/output/mixer/

mkdir -p /exports/eddie/scratch/s2093301/MDD_MCP/output/{err,out,mixer}

qsub -V -e ./output/err/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID_\$TASK_ID -o ./output/out/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID_\$TASK_ID $SCRIPT_DIR/bivariate_mixer_master.sh
qsub -V - ./output/err/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID_\$TASK_ID -o ./output/out/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID_\$TASK_ID $SCRIPT_DIR/combine_plot__master.sh

# MCP X MDD (EUR EHR)
export RUN_UNIVARIATE_1=0 #MCP already run

export SUMSTATS_FILE2=MDD_EHR_MiXeR.sumstats.gz
export TRAIT2=MDD_EHR

qsub -V -e ./output/err/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID_\$TASK_ID -o ./output/out/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID_\$TASK_ID $SCRIPT_DIR/bivariate_mixer_master.sh
qsub -V - ./output/err/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID_\$TASK_ID -o ./output/out/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID_\$TASK_ID $SCRIPT_DIR/combine_plot__master.sh


# MCP X MDD (EUR Quest)
export SUMSTATS_FILE2=MDD_quest_MiXeR.sumstats.gz
export TRAIT2=MDD_quest

qsub -V -e ./output/err/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID_\$TASK_ID -o ./output/out/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID_\$TASK_ID $SCRIPT_DIR/bivariate_mixer_master.sh
qsub -V - ./output/err/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID_\$TASK_ID -o ./output/out/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID_\$TASK_ID $SCRIPT_DIR/combine_plot__master.sh
