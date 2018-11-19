rule filtlong:
    input:
        "data/{run}/porechopped/{sample}.fastq.gz"
    output:
        "data/{run}/filtlong/{sample}_filtlong.fastq.gz"
    threads:
        cluster_config["filtlong"]["nCPUs"]
    resources:
        mem_mb=cluster_config["filtlong"]["memory"]
    params:
        min_read_length = config["min_read_length"],
        keep_percent = config["keep_percent"],
        mean_q_weight = config["mean_q_weight"]
    log:
        "logs/filtlong_{run}_{sample}.log"
    singularity:
        config["container"]
    shell:
        "filtlong  "
