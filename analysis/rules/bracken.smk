rule build_bracken_db:
    input:
        rules.build_kraken2_db.output
    output:
        "data/kraken2_db/database1000mers.kmer_distrib"
    threads:
        cluster_config["build_bracken_db"]["nCPUs"]
    resources:
        mem_mb = cluster_config["build_bracken_db"]["memory"]
    params:
        kraken2_db = "data/kraken2_db",
        read_length = config["min_read_length"]
    singularity:
        config["container"]
    log:
        "logs/build_bracken_db.log"
    shell:
        """
        bracken-build -d {params.kraken2_db} \
          -threads {threads} \
          -l {params.read_length} 2> {log}
        """

rule bracken_classify:
    input:
        rules.build_bracken_db.output,
        report = "data/{run}/kraken2/kraken2_classification_{run}_{sample}.kreport"
    output:
        "data/{run}/bracken/bracken_classification_{run}_{sample}.bracken"
    resources:
        mem_mb = cluster_config["bracken_classify"]["memory"]
    singularity:
        config["container"]
    params:
        db = "data/kraken2_db",
        read_length = config["min_read_length"],
        threshold = config["bracken_threshold"]
    log:
        "logs/bracken_classify_{run}_{sample}.log"
    shell:
        """
        bracken -d {params.db} \
          -i {input.report} \
          -o {output} \
          -r {params.read_length} \
          -t {params.threshold}
        """
