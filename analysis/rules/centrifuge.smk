rule download_centrifuge_db:
    output:
        expand("data/centrifuge_db/p_compressed.{n}.cf", n=[1, 2, 3, 4]),
    threads: 1
    resources:
        mem_mb = 500
    params:
        db_url = config["download_centrifuge_db"]["url"],
        output = "data/centrifuge_db/" + config["download_centrifuge_db"]["url"].split("/")[-1],
        md5_hash = config["download_centrifuge_db"]["md5"]
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
        fastq = "data/{run}/filtlong/{sample}.filtered.fastq.gz"
    output:
        report = "data/{run}/centrifuge/centrifuge_report_{run}_{sample}.tsv",
        classification = "data/{run}/centrifuge/centrifuge_classification_{run}_{sample}.tab"
    threads:
        config["centrifuge"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["centrifuge"]["memory"]
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
    threads: 1
    resources:
        mem_mb = 500
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
        expand("data/centrifuge_16s_db/silva_16s.{db_idx}.cf", db_idx=range(1, 5))
    threads:
        config["build_centrifuge_16s_db"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["build_centrifuge_16s_db"]["memory"]
    params:
        prefix = "data/centrifuge_16s_db/silva_16s"
    log:
        "logs/build_centrifuge_16s_db.log"
    singularity:
        config["container"]
    shell:
        """
        centrifuge-build \
          --threads {threads} \
          --conversion-table {input.conversion_table} \
          --taxonomy-tree {input.tax_tree} \
          --name-table {input.name_table} \
          {input.ref_seqs} \
          {params.prefix} 2> {log}
        """

rule centrifuge_16s_classify:
    input:
        rules.build_centrifuge_16s_db.output,
        fastq = "data/{run}/filtlong/{sample}.filtered.fastq.gz"
    output:
        report = "data/{run}/centrifuge/centrifuge_16s_report_{run}_{sample}.tsv",
        classification = "data/{run}/centrifuge/centrifuge_16s_classification_{run}_{sample}.tab"
    threads:
        config["centrifuge_16s"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["centrifuge_16s"]["memory"]
    params:
        index_prefix = "data/centrifuge_16s_db/silva_16s"
    singularity:
        config["container"]
    log:
        "logs/centrifuge_16s_classify_{run}_{sample}.log"
    shell:
        """
        centrifuge -x {params.index_prefix} \
          -U {input.fastq} \
          --threads {threads} \
          --report-file {output.report} \
          -S {output.classification} \
          --met-stderr 2> {log}
        """

rule centrifuge_16s_krakenstyle_report:
    input:
        "data/{run}/centrifuge/centrifuge_16s_classification_{run}_{sample}.tab"
    output:
        "data/{run}/centrifuge/centrifuge_16s_classification_{run}_{sample}.kreport"
    threads: 1
    resources:
        mem_mb = 500
    params:
        index_prefix = "data/centrifuge_16s_db/silva_16s"
    log:
        "logs/centrifuge_16s_kreport_{run}_{sample}.log"
    singularity:
        config["container"]
    shell:
        """
        centrifuge-kreport -x {params.index_prefix} {input} > {output} 2> {log}
        """

rule build_centrifuge_full_kraken_db:
    input:
        rules.build_kraken2_db.output
    output:
        expand("data/centrifuge_full_kraken_db/kraken_full.{db_idx}.cf", db_idx=range(1, 5))
    threads:
        config["build_centrifuge_full_kraken_db"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["build_centrifuge_full_kraken_db"]["memory"]
    params:
        prefix = "data/centrifuge_full_kraken_db/kraken_full",
        conversion_table = "data/kraken2_db/seqid2taxid.map",
        tax_tree = "data/kraken2_db/taxonomy/nodes.dmp",
        name_table = "data/kraken2_db/taxonomy/names.dmp",
        ref_seqs = "data/kraken2_db/library/bacteria/library.fna"
    log:
        "logs/build_centrifuge_full_kraken_db.log"
    singularity:
        config["container"]
    shell:
        """
        centrifuge-build \
          --threads {threads} \
          --conversion-table {params.conversion_table} \
          --taxonomy-tree {params.tax_tree} \
          --name-table {params.name_table} \
          {params.ref_seqs} \
          {params.prefix} 2> {log}
        """

rule centrifuge_full_kraken_classify:
    input:
        rules.build_centrifuge_full_kraken_db.output,
        fastq = "data/{run}/filtlong/{sample}.filtered.fastq.gz"
    output:
        report = "data/{run}/centrifuge/centrifuge_full_kraken_report_{run}_{sample}.tsv",
        classification = "data/{run}/centrifuge/centrifuge_full_kraken_classification_{run}_{sample}.tab"
    threads:
        config["centrifuge"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["centrifuge"]["memory"]
    params:
        index_prefix = "data/centrifuge_full_kraken_db/kraken_full"
    singularity:
        config["container"]
    log:
        "logs/centrifuge_full_kraken_classify_{run}_{sample}.log"
    shell:
        """
        centrifuge -x {params.index_prefix} \
          -U {input.fastq} \
          --threads {threads} \
          --report-file {output.report} \
          -S {output.classification} \
          --met-stderr 2> {log}
        """

rule centrifuge_full_kraken_krakenstyle_report:
    input:
        "data/{run}/centrifuge/centrifuge_full_kraken_classification_{run}_{sample}.tab"
    output:
        "data/{run}/centrifuge/centrifuge_full_kraken_classification_{run}_{sample}.kreport"
    threads: 1
    resources:
        mem_mb = 500
    params:
        index_prefix = "data/centrifuge_full_kraken_db/kraken_full"
    log:
        "logs/centrifuge_full_kraken_krakenstyle_report_{run}_{sample}.log"
    singularity:
        config["container"]
    shell:
        """
        centrifuge-kreport -x {params.index_prefix} {input} > {output} 2> {log}
        """
