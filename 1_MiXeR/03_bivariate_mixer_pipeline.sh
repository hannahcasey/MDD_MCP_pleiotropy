export APPTAINER_BIND=/exports/eddie/scratch/s2093301/MDD_MCP/
export MIXER_DIR=/exports/eddie/scratch/s2093301/MDD_MCP/resources/mixer/
export SCRIPT_DIR=/exports/eddie/scratch/s2093301/MDD_MCP/bin
export SUMSTATS_DIR=/exports/eddie/scratch/s2093301/MDD_MCP/resources/sumstats/

export RESULTS_DIR=/exports/eddie/scratch/s2093301/MDD_MCP/output/mixer/
mkdir -p /exports/eddie/scratch/s2093301/MDD_MCP/output/{err,out,mixer}

# MCP X MDD (EUR meta-analysis - UKB and 23andMe removed)
export SUMSTATS_FILE1=MCP_MiXeR.sumstats.gz
export TRAIT1=MCP
export RUN_UNIVARIATE_1=1

export SUMSTATS_FILE2=MDD_MiXeR.sumstats.gz
export TRAIT2=MDD
export RUN_UNIVARIATE_2=1

JOB1=$(qsub -V \
    -e ./output/err/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID \
    -o ./output/out/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID \
    $SCRIPT_DIR/03_bivariate_mixer_master.sh)

PLOT1=$(qsub -V -W depend=afterok:${JOB1} \
    -e ./output/err/bivariate_mixer_plot_${TRAIT1}_${TRAIT2}_\$JOB_ID \
    -o ./output/out/bivariate_mixer_plot_${TRAIT1}_${TRAIT2}_\$JOB_ID \
    $SCRIPT_DIR/03_combine_plot_master.sh)

# MCP X MDD (EUR EHR)
export RUN_UNIVARIATE_1=0 #MCP already run

export SUMSTATS_FILE2=MDD_EHR_MiXeR.sumstats.gz
export TRAIT2=MDD_EHR

JOB2=$(qsub -V -W depend=afterok:${PLOT1} \
    -e ./output/err/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID \
    -o ./output/out/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID \
    $SCRIPT_DIR/03_bivariate_mixer_master.sh)

PLOT2=$(qsub -V -W depend=afterok:${JOB2} \
    -e ./output/err/bivariate_mixer_plot_${TRAIT1}_${TRAIT2}_\$JOB_ID \
    -o ./output/out/bivariate_mixer_plot_${TRAIT1}_${TRAIT2}_\$JOB_ID \
    $SCRIPT_DIR/03_combine_plot_master.sh)

# MCP X MDD (EUR Quest)
export SUMSTATS_FILE2=MDD_quest_MiXeR.sumstats.gz
export TRAIT2=MDD_quest

JOB3=$(qsub -V -W depend=afterok:${PLOT2} \
    -e ./output/err/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID \
    -o ./output/out/bivariate_mixer_${TRAIT1}_${TRAIT2}_\$JOB_ID \
    $SCRIPT_DIR/03_bivariate_mixer_master.sh)

PLOT3=$(qsub -V -W depend=afterok:${JOB3} \
    -e ./output/err/bivariate_mixer_plot_${TRAIT1}_${TRAIT2}_\$JOB_ID \
    -o ./output/out/bivariate_mixer_plot_${TRAIT1}_${TRAIT2}_\$JOB_ID \
    $SCRIPT_DIR/03_combine_plot_master.sh)
