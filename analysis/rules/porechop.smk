if MULTIPLEXED:
    INPUT = "data/basecalled/"
else:
    INPUT = expand("data/basecalled/{sample}.fastq.gz", sample=SAMPLES)


def determine_output_format(wildcards, output):
    if MULTIPLEXED:
        result = "--barcode_dir data/porechopped"
    else:
        result = "--output {output}".format(output=output[0])
    return result


rule porechop:
    input:
        INPUT
    output:
        expand("data/porechopped/{sample}.fastq.gz", sample=SAMPLES)
    threads:
        cluster_config["porechop"]["nCPUs"]
    resources:
        mem_mb=cluster_config["porechop"]["memory"]
    params:
        output_type=determine_output_format,
        unassigned=("--discard_unassigned" if MULTIPLEXED else "")
    log:
        "logs/porechop.log"
    singularity:
        config["container"]
    shell:
        "porechop --input {input}  {params.output_type} --threads {threads} "
        "--check_reads 25000 --extra_end_trim 10 --discard_middle "
        "{params.unassigned} --format fastq.gz > {log}"
