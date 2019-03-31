#!/usr/bin/env bash
CLUSTER_CMD=("bsub -n {threads} -R \"select[mem>{resources.mem_mb}] rusage[mem={resources.mem_mb}] span[hosts=1]\" -M {resources.mem_mb} -o {cluster.output} -e {cluster.error} -J {cluster.name} -q {cluster.queue}")
JOB_NAME="$1"
LOG_DIR="logs/"

bsub -R "select[mem>1500] rusage[mem=1500]" \
    -M 1500 \
    -o "$LOG_DIR"/cluster_"$JOB_NAME".o \
    -e "$LOG_DIR"/cluster_"$JOB_NAME".e \
    -J "$JOB_NAME" \
    -q research-rh74 \
    snakemake --use-singularity \
    --cluster-config cluster.yaml \
    --jobs 2000 \
    --restart-times 3 \
    --cluster "${CLUSTER_CMD[@]}"

exit 0
