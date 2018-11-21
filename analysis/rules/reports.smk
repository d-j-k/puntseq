rule plot_pre_filtering:
    input:
        "data/porechopped/{sample}.fastq.gz"
    output:
        "data/plots/{sample}_pre_filtering.pdf"
    log:
        "logs/pistis_pre_filtering_{sample}.log"
    resources:
        mem_mb=cluster_config["plot_pre_filtering"]["memory"]
    params:
        downsample='--downsample 0'
    singularity:
        config["container"]
    shell:
        "pistis --fastq {input} --output {output} {params.downsample} 2> {log} "


rule plot_post_filtering:
    input:
        fastq="data/filtered/{sample}_filtered.fastq.gz",
        bam="data/sorted/{sample}_sorted.bam"
    output:
        "data/plots/{sample}_post_filtering.pdf"
    log:
        "logs/pistis_post_filtering_{sample}.log"
    resources:
        mem_mb=cluster_config["plot_post_filtering"]["memory"]
    params:
        downsample='--downsample 0'
    singularity:
        config["container"]
    shell:
        "pistis --fastq {input.fastq} --output {output} --bam {input.bam} "
        "{params.downsample} 2> {log} "


rule stats_pre_filtering:
    input:
        "data/porechopped/{sample}.fastq.gz"
    output:
        "data/stats/{sample}_pre_filtering.txt"
    log:
        "logs/nanostat_pre_filtering_{sample}.log"
    threads:
        cluster_config["stats_pre_filtering"]["nCPUs"]
    resources:
        mem_mb=cluster_config["stats_pre_filtering"]["memory"]
    singularity:
        config["container"]
    shell:
        "NanoStat --fastq {input} --name {output} --threads {threads} "
        "--readtype 1D 2> {log}"

rule stats_post_filtering:
    input:
        "data/filtered/{sample}_filtered.fastq.gz"
    output:
        "data/stats/{sample}_post_filtering.txt"
    log:
        "logs/nanostat_post_filtering_{sample}.log"
    threads:
        cluster_config["stats_post_filtering"]["nCPUs"]
    resources:
        mem_mb=cluster_config["stats_post_filtering"]["memory"]
    singularity:
        config["container"]
    shell:
        "NanoStat --fastq {input} --name {output} --threads {threads} "
        "--readtype 1D 2> {log}"


rule report:
    input:
        plot_pre_filter="data/plots/{sample}_pre_filtering.pdf",
        plot_post_filter="data/plots/{sample}_post_filtering.pdf",
        stats_pre_filter="data/stats/{sample}_pre_filtering.txt",
        stats_post_filter="data/stats/{sample}_post_filtering.txt",
        mykrobe="data/mykrobe/{sample}/{sample}_predict.json",
        porechop_log="logs/porechop.log"
    output:
        "docs/report_{sample}.html"
    log:
        "logs/report_{sample}.log"
    params:
        sample="{sample}"
    script:
        "../scripts/report.py"
