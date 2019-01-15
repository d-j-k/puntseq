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
        "data/{run}/centrifuge/centrifuge_classification_{run}_{sample}.kreport"
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

rule build_centrifuge_16s_db:
    input:
        name_table = rules.build_kraken2_16s_db.output.name_table,
        tax_tree = rules.build_kraken2_16s_db.output.tax_tree,
        conversion_table = rules.build_kraken2_16s_db.output.conversion_table,
        ref_seqs = rules.build_kraken2_16s_db.output.ref_seqs
    output:
        "data/centrifuge_16s_db/silva_16s.1.cf",
        "data/centrifuge_16s_db/silva_16s.2.cf",
        "data/centrifuge_16s_db/silva_16s.3.cf"
    threads:
        cluster_config["build_centrifuge_16s_db"]["nCPUs"]
    resources:
        mem_mb = cluster_config["build_centrifuge_16s_db"]["memory"]
    params:
        prefix = "data/centrifuge_16s_db/silva_16s"
    log:
        "logs/build_centrifuge_16s_db.log"
    singularity:
        config["container"]
    shell:
        """
        centrifuge-build -c {input.ref_seqs} \
          --threads {threads} \
          --conversion-table {input.conversion_table} \
          --taxonomy-tree {input.tax_tree} \
          --name-table {input.name_table} \ 
          {params.prefix} 2> {log}
        """
