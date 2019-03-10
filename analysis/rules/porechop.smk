def determine_demultiplex_action(wildcards, input, output, threads, resources):
    df = samples.loc[(wildcards.run)]
    out_dir = Path(
        "data/{run}/porechopped/".format(run=wildcards.run)
    )

    expected_barcodes = df["barcode"]
    is_multiplexed = not any(expected_barcodes.isnull())

    if is_multiplexed:
        out = out_dir
    else:
        out = out_dir / "{}.{}".format(wildcards.run, config["porechop_out_format"])

    result = {
        "df": df,
        "out_dir": out_dir,
        "out": out,
        "is_multiplexed": is_multiplexed,
        "expected_barcodes": set(expected_barcodes),
    }
    return result



rule porechop:
    input:
        fastq  = "data/{run}/basecalled/{run}_all_passed.fastq.gz"
    output:
        "data/{run}/porechopped/DEMULTIPLEX_COMPLETE"
    threads:
        config["porechop"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["porechop"]["memory"]
    singularity:
        config["container"]
    params:
        option = determine_demultiplex_action,
        check_reads = config["porechop"]["check_reads"],
        out_format = config["porechop"]["output_format"],
    log:
        "logs/porechop_{run}.log"
    shell:
        """
        bash analysis/scripts/porechop.sh {params.option[is_multiplexed]} \
            {input.fastq} \
            {params.option[out]}
            {params.option[classification_path]} \
            {output} \
            {threads} \
            {params.check_reads} \
            {params.out_format} 2> {log}
        """


rule fix_filenames:
    input:
        "data/{run}/porechopped/DEMULTIPLEX_COMPLETE"
    output:
        "data/{run}/porechopped/{sample}.trimmed.fastq.gz"
    threads: 1
    resources:
        mem_mb = 250
    params:
        option = determine_demultiplex_action
    log:
        "logs/fix_filenames_{run}.log"
    run:
        original_path = Path(output[0].split(".")[0] + config["porechop"]["output_format"])
        if not original_path.is_file():
            raise FileNotFoundError("Expected porechop output {} not found.".format(original_path))

        original_path.replace(output[0])
