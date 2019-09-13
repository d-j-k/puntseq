#rule filtlong:
#    input:
#        "data/{run}/porechopped/{sample}.trimmed.fastq.gz"
#    output:
#        "data/{run}/filtlong/{sample}.filtered.fastq.gz"
#    threads: 1
#    resources:
#        mem_mb = lambda wildcards, attempt: attempt * config["filtlong"]["memory"]
#    params:
#        min_read_length = config["filtlong"]["min_len"],
#    log:
#        "logs/filtlong_{run}_{sample}.log"
#    singularity:
#        config["container"]
#    shell:
#        """
#        filtlong --min_length {params.min_read_length} \
#            --verbose {input} 2> {log} | gzip > {output}
#        """

rule filter_length:
    input:
        "data/{run}/porechopped/{sample}.trimmed.fastq.gz"
    output:
        "data/{run}/filtlong/{sample}.filtered.fastq.gz"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["filter_length"]["memory"]
    params:
        min_len = config["filter_length"]["min_len"],
        max_len = config["filter_length"]["max_len"],
    log:
        "logs/filter_length_{run}_{sample}.log"
    singularity: 
        "containers/nanofilt.v2.5.0.simg"
    shell:
        """
        gzip -d -c {input} | \
        NanoFilt --length {params.min_len} --maxlength {params.max_len} | \
        gzip > {output} 2> {log}
        """
