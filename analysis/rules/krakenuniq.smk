rule download_krakenuniq_db:
    output:
        "data/krakenuniq_db/taxonomy/names.dmp",
        "data/krakenuniq_db/taxonomy/nodes.dmp"
    threads:
        cluster_config["download_krakenuniq_db"]["nCPUs"]
    resources:
        mem_mb = cluster_config["download_krakenuniq_db"]["memory"]
    params:
        taxa = "archaea,bacteria,viral,fungi,protozoa",
        db_dir = "data/krakenuniq_db",
        database = config["krakenuniq_database"]
    singularity:
        config["container"]
    log:
        "logs/krakenuniq_db_download.log"
    shell:
        """
        krakenuniq-download --db {params.db_dir} \
           --taxa {params.taxa} \
          --threads {threads} \
          --dust \
          {params.database} 2> {log}

        """

rule krakenuniq:
    input:
        rules.download_krakenuniq_db.output,
        fastq = "data/{run}/filtlong/{sample}_filtered.fastq.gz"
    output:
        report = "data/{run}/krakenuniq/krakenuniq_classification_{run}_{sample}.kreport",
        outfile = "data/{run}/krakenuniq/krakenuniq_classification_{run}_{sample}.out"
    threads:
        cluster_config["krakenuniq"]["nCPUs"]
    resources:
        mem_mb = cluster_config["krakenuniq"]["memory"]
    params:
        db_dir = "data/krakenuniq_db"
    shell:
        """
        krakenuniq --db {params.db_dir} \
          --report-file {output.report} \
          --output {output.outfile} \
          --threads {threads}
        """
