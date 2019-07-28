rule plot_post_filtering:
    input:
        "data/{run}/filtlong/{sample}.filtered.fastq.gz"
    output:
        "data/{run}/plots/{sample}.filtered.pdf"
    log:
        "logs/pistis_{run}_{sample}.log"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["plot"]["memory"]
    params:
        downsample = config["plot"]["downsample"]
    shell:
        """
        pistis --fastq {input} --output {output} \
          --downsample {params.downsample} 2> {log}
        """


rule stats_post_filtering:
    input:
        "data/{run}/filtlong/{sample}.filtered.fastq.gz"
    output:
        "data/{run}/stats/{sample}.stats.txt"
    log:
        "logs/nanostat_{run}_{sample}.log"
    threads:
        config["stats"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["stats"]["memory"]
    shell:
        """
        NanoStat --fastq {input} --name {output} --threads {threads} 2> {log}
        """
