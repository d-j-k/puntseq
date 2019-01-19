bracken_db_output = "data/kraken2_db/database{}mers.kmer_distrib".format(config["min_read_length"])

rule build_bracken_db:
    input:
        rules.build_kraken2_db.output
    output:
        bracken_db_output
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
          -t {threads} \
          -l {params.read_length} 2> {log}
        """

bracken_16s_db_output = "data/kraken2_16s_db/database{}mers.kmer_distrib".format(config["min_read_length"])

rule build_bracken_16s_db:
    input:
        rules.build_kraken2_16s_db.output
    output:
        bracken_16s_db_output
    params:
        kraken2_db = "data/kraken2_16s_db",
        read_length = config["min_read_length"]
    threads:
        cluster_config["build_bracken_16s_db"]["nCPUs"]
    resources:
        mem_mb = cluster_config["build_bracken_16s_db"]["memory"]
    log:
        "logs/build_bracken_16s_db.log"
    singularity:
        config["container"]
    shell:
        """
        bracken-build -d {params.kraken2_db} \
          -t {threads} \
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
          -t {params.threshold} 2> {log}
        """

rule bracken_16s_classify:
    input:
        rules.build_bracken_16s_db.output,
        report = "data/{run}/kraken2/kraken2_16s_classification_{run}_{sample}.kreport"
    output:
        "data/{run}/bracken/bracken_16s_classification_{run}_{sample}.bracken"
    resources:
        mem_mb = cluster_config["bracken_16s_classify"]["memory"]
    singularity:
        config["container"]
    params:
        db = "data/kraken2_16s_db",
        read_length = config["min_read_length"],
        threshold = config["bracken_threshold"]
    log:
        "logs/bracken_16s_classify_{run}_{sample}.log"
    shell:
        """
        bracken -d {params.db} \
          -i {input.report} \
          -o {output} \
          -r {params.read_length} \
          -t {params.threshold} 2> {log}
        """
