bracken_db_output = "data/kraken2_db/database{}mers.kmer_distrib".format(config["min_read_length"])

rule build_bracken_db:
    input:
        rules.build_kraken2_db.output
    output:
        bracken_db_output
    threads:
        config["build_bracken_db"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["build_bracken_db"]["memory"]
    params:
        kraken2_db = "data/kraken2_db",
        read_length = config["build_bracken_db"]["read_length"],
        kmer_length = config["build_bracken_db"]["kmer_length"],
    singularity:
        config["container"]
    log:
        "logs/build_bracken_db.log"
    shell:
        """
        bracken-build -d {params.kraken2_db} \
          -t {threads} \
          -k {params.kmer_length} \
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
        read_length = config["build_bracken_16s_db"]["read_length"],
        kmer_length = config["build_bracken_16s_db"]["kmer_length"],
    threads:
        config["build_bracken_16s_db"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["build_bracken_16s_db"]["memory"]
    log:
        "logs/build_bracken_16s_db.log"
    singularity:
        config["container"]
    shell:
        """
        bracken-build -d {params.kraken2_db} \
          -t {threads} \
          -k {params.kmer_length} \
          -l {params.read_length} 2> {log}
        """

rule bracken_classify:
    input:
        rules.build_bracken_db.output,
        report = "data/{run}/kraken2/kraken2_classification_{run}_{sample}.kreport"
    output:
        "data/{run}/bracken/bracken_classification_{run}_{sample}.bracken"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["bracken"]["memory"]
    singularity:
        config["container"]
    params:
        db = "data/kraken2_db",
        read_length = config["bracken"]["read_length"],
        threshold = config["bracken"]["threshold"]
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
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["bracken"]["memory"]
    singularity:
        config["container"]
    params:
        db = "data/kraken2_16s_db",
        read_length = config["bracken"]["read_length"],
        threshold = config["bracken"]["threshold"]
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

bracken_16s_db_k21_output = "data/kraken2_16s_db_k21/database{}mers.kmer_distrib".format(config["min_read_length"])

rule build_bracken_16s_db_k21:
    input:
        rules.build_kraken2_16s_db_k21.output
    output:
        bracken_16s_db_k21_output
    params:
        kraken2_db = "data/kraken2_16s_db_k21",
        read_length = config["build_bracken_16s_db_k21"]["read_length"],
        kmer_length = config["build_bracken_16s_db_k21"]["kmer_length"],
    threads:
        config["build_bracken_16s_db_k21"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["build_bracken_16s_db_k21"]["memory"]
    log:
        "logs/build_bracken_16s_db_k21.log"
    singularity:
        config["container"]
    shell:
        """
        bracken-build -d {params.kraken2_db} \
          -t {threads} \
          -k {params.kmer_length} \
          -l {params.read_length} 2> {log}
        """

rule bracken_16s_k21_classify:
    input:
        rules.build_bracken_16s_db_k21.output,
        report = "data/{run}/kraken2/kraken2_16s_k21_classification_{run}_{sample}.kreport"
    output:
        "data/{run}/bracken/bracken_16s_k21_classification_{run}_{sample}.bracken"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["bracken"]["memory"]
    singularity:
        config["container"]
    params:
        db = "data/kraken2_16s_db_k21",
        read_length = config["bracken"]["read_length"],
        threshold = config["bracken"]["threshold"]
    log:
        "logs/bracken_16s_k21_classify_{run}_{sample}.log"
    shell:
        """
        bracken -d {params.db} \
          -i {input.report} \
          -o {output} \
          -r {params.read_length} \
          -t {params.threshold} 2> {log}
        """
