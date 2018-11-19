
def determine_output_format(wildcards, output):
    if MULTIPLEXED:
        result = "--barcode_dir data/porechopped"
    else:
        result = "--output {output}".format(output=output[0])
    return result


rule porechop:
    input:
        expand("data/basecalled/{sample}.fastq.gz", sample=SAMPLES)
    output:
        expand("data/porechopped/{sample}.fastq.gz", sample=SAMPLES)
    threads:
        cluster_config["porechop"]["nCPUs"]
    resources:
        mem_mb=cluster_config["porechop"]["memory"]
    params:
        barcode_dir = "data/{wildcars.run}/porechopped",
        check_reads = config["check_reads"]
    log:
        "logs/porechop.log"
    singularity:
        config["container"]
    shell:
        "porechop --input {input}  {params.barcode_dir} --threads {threads} "
        "--check_reads {params.check_reads} --extra_end_trim 10 --discard_middle "
        "--discard_unassigned --format fastq.gz > {log}"
