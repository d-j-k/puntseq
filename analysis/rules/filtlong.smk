rule filtlong:
    input:
        "data/{run}/porechopped/{sample}.trimmed.fastq.gz"
    output:
        "data/{run}/filtlong/{sample}.filtered.fastq.gz"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["filtlong"]["memory"]
    params:
        min_read_length = config["filtlong"]["min_len"],
    log:
        "logs/filtlong_{run}_{sample}.log"
    singularity:
        config["container"]
    shell:
        """
        filtlong --min_length {params.min_read_length} \
            --verbose {input} 2> {log} | gzip > {output}
        """
