rule map:
    input:
        target = rules.convert_silva_db_to_dna.output,
        query = "data/{run}/filtlong/{sample}.filtered.fastq.gz"
    output:
        "data/{run}/map/{sample}.filtered.bam"
    threads: 
        config["map"].get("threads", 1)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["map"]["memory"]
    #singularity:
    #    config["container"]
    params:
        extra = "-K 25M --no-kalloc --print-qname -aLx map-ont"
    log:
        "logs/map/{run}_{sample}.log"
    shell:
        """
        minimap2 {params.extra} -t {threads} {input.target} {input.query} | \
            samtools sort -@{threads} -o {output} - 2> {log}
        """
