rule plot_post_filtering:
    input:
        "data/{run}/filtlong/{sample}_filtered.fastq.gz"
    output:
        "data/{run}/plots/{sample}_filtered.pdf"
    log:
        "logs/pistis_{run}_{sample}.log"
    resources:
        mem_mb = cluster_config["plot_post_filtering"]["memory"] 
    params:
        downsample = config["plot_downsample"]
    shell:
        """
        pistis --fastq {input} --output {output} \
          --downsample {params.downsample} 2> {log}
        """


rule stats_post_filtering:
    input:
        "data/{run}/filtlong/{sample}_filtered.fastq.gz"
    output:
        "data/{run}/stats/{sample}_stats.txt"
    log:
        "logs/nanostat_{run}_{sample}.log"
    threads: 
        cluster_config["stats_post_filtering"]["nCPUs"]
    resources:
        mem_mb = cluster_config["stats_post_filtering"]["memory"]
    shell:
        """
        NanoStat --fastq {input} --name {output} --threads {threads} 2> {log}
        """

