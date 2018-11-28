rule download_kraken2_taxonomy:
    output:
        "data/kraken2_db/taxonomy/names.dmp",
        "data/kraken2_db/taxonomy/nodes.dmp"
    singularity:
        config["container"]
    params:
        db = "data/kraken2_db"
    log:
        "logs/download_kraken2_taxonomy.log"
    shell:
        """
        kraken2-build --download-taxonomy --db {params.db} 2> {log}
        """

rule download_kraken2_db:
    input:
        rules.download_kraken2_taxonomy.output
    output:
        "data/kraken2_db/library/bacteria/library.fna",
        "data/kraken2_db/library/archaea/library.fna"
    threads:
        cluster_config["build_kraken2_db"]["nCPUs"]
    resources:
        mem_mb = cluster_config["build_kraken2_db"]["memory"]
    singularity:
        config["container"]
    log:
        "logs/download_kraken2_db.log"
    params:
        db = "data/kraken2_db"
    shell:
        """
        kraken2-build --download-library archaea --db {params.db} --threads {threads} 2>> {log}
        kraken2-build --download-library bacteria --db {params.db} --threads {threads} 2>> {log}
        """

rule build_kraken2_db:
    input:
        rules.download_kraken2_db.output
    output:
        "data/kraken2_db/hash.k2d",
        "data/kraken2_db/opts.k2d",
        "data/kraken2_db/taxo.k2d"
    resources:
        mem_mb = cluster_config["build_kraken2_db"]["memory"]
    threads:
        cluster_config["build_kraken2_db"]["nCPUs"]
    singularity:
        config["container"]
    params:
        db = "data/kraken2_db"
    log:
        "logs/build_kraken2_db.log"
    shell:
        """
        kraken2-build --build --db {params.db} 2> {log}
        touch deleteme
        """

rule kraken2_classify:
    input:
        rules.build_kraken2_db.output,
        fastq = "data/{run}/filtlong/{sample}_filtered.fastq.gz"
    output:
        report = "data/{run}/kraken2/kraken2_classification_{run}_{sample}.kreport",
        outfile = "data/{run}/kraken2/kraken2_classification_{run}_{sample}.out"
    threads:
        cluster_config["kraken2_classify"]["nCPUs"]
    resources:
        mem_mb = cluster_config["kraken2_classify"]["memory"]
    params:
        db_dir = "data/kraken2_db"
    log:
        "logs/kraken2_classify_{run}_{sample}.log"
    singularity:
        config["container"]
    shell:
        """
        kraken2 --db {params.db_dir} \
          --output {output.outfile} \
          --report {output.report} \
          --fastq-input \
          --gzip-compressed \
          --threads {threads} {input.fastq} 2> {log}
        """
