rule download_centrifuge_db:
    output:
        expand("data/centrifuge_db/p_compressed.{n}.cf", n=[1, 2, 3, 4]),
    params:
        db_url = config["centrifuge_db_url"],
        output = "data/centrifuge_db/" + config["centrifuge_db_url"].split("/")[-1],
        md5_hash = config["db_md5"]
    log:
        "logs/download_centrifuge_db.log"
    shell:
        """
        scripts/download_centrifuge_db.sh \
          {params.db_url} \
          {params.output} \
          {params.md5_hash} 2> {log}
        """

rule centrifuge:
    input:
        rules.download_centrifuge_db.output,
        fastq = "data/{run}/filtlong/{sample}_filtered.fastq.gz"
    output:
        report = "data/{run}/centrifuge/centrifuge_report_{run}_{sample}.tsv",
        classification = "data/{run}/centrifuge/centrifuge_classification_{run}_{sample}.tab"
    threads:
        cluster_config["centrifuge"]["nCPUs"]
    resources:
        mem_mb = cluster_config["centrifuge"]["memory"]
    params:
        index_prefix = "data/centrifuge_db/p_compressed"
    singularity:
        config["container"]
    log:
        "logs/centrifuge_{run}_{sample}.log"
    shell:
        """
        centrifuge -x {params.index_prefix} \
          -U {input.fastq} \
          --threads {threads} \
          --report-file {output.report} \
          -S {output.classification} \
          --met-stderr 2> {log}
        """

rule centrifuge_krakenstyle_report:
    input:
        "data/{run}/centrifuge/centrifuge_classification_{run}_{sample}.tab"
    output:
        "data/{run}/centrifuge/centrifuge_classification_kreport_{run}_{sample}.tab"
    params:
        index_prefix = "data/centrifuge_db/p_compressed"
    log:
        "logs/centrifuge_kreport_{run}_{sample}.log"
    singularity:
        config["container"]
    shell:
        """
        centrifuge-kreport -x {params.index_prefix} {input} > {output} 2> {log}
        """
