rule filtlong:
    input:
        rules.porechop.output
    output:
        "data/{run}/filtlong/{sample}_filtered.fastq.gz"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000
    params:
        min_read_length = config["min_read_length"],
        keep_percent = config["keep_percent"],
        mean_q_weight = config["mean_q_weight"]
    log:
        "logs/filtlong_{run}_{sample}.log"
    singularity:
        config["container"]
    shell:
        """
        filtlong --min_length {params.min_read_length} \
            --keep_percent {params.keep_percent} \
            --mean_q_weight {params.mean_q_weight} \
            --verbose {input} 2> {log} | gzip > {output}
        """
