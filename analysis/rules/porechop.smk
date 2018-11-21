rule porechop:
    input:
        "data/{run}/basecalled/{run}_all_passed.fastq.gz"
    output:
        expand("data/{{run}}/porechopped/{sample}.fastq.gz", sample=SAMPLES)
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000
    params:
        barcode_dir = "data/{wildcards.run}/porechopped",
        check_reads = config["check_reads"],
        extra_end_trim = config["extra_end_trim"],
        out_format = config["porechop_out_format"]
    log:
        "logs/porechop_{run}.log"
    singularity:
        config["container"]
    shell:
        """
        porechop --input {input} \
          --barcode_dir {params.barcode_dir} \
          --threads {threads} \
          --check_reads {params.check_reads} \
          --extra_end_trim {params.extra_end_trim} \
          --discard_middle \
          --discard_unassigned \
          --format {params.out_format} &> {log}
        """
