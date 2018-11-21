#!/usr/bin/env bash
CLUSTER_CMD=("bsub -n {cluster.nCPUs} -R {cluster.resources} -M {cluster.memory} -o {cluster.output} -e {cluster.error} -J {cluster.name}")
JOB_NAME="$1"

bsub -R "select[mem>1000] rusage[mem=1000]" \
  -M 1000 \
  -o logs/cluster_"$JOB_NAME".o \
  -e logs/cluster_"$JOB_NAME".e \
  -J "$JOB_NAME" \
  snakemake --use-singularity \
    --cluster-config cluster.yaml \
    --jobs 500 \
    --cluster "${CLUSTER_CMD[@]}" \
    --restart-times 3

exit 0
