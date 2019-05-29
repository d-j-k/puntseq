rule download_centrifuge_taxonomy:
    output:
        name_table = "data/centrifuge_db/taxonomy/names.dmp",
        tax_tree = "data/centrifuge_db/taxonomy/nodes.dmp",
    threads: 1
    resources:
        mem_mb = 200
    singularity:
        config["container"]
    log:
        "logs/download_centrifuge_taxonomy.log"
    shell:
        """
        centrifuge-download -o data/centrifuge_db/taxonomy taxonomy 2> {log}
        """


rule download_centrifuge_library:
    output:
        mapping = "data/centrifuge_db/seqid2taxid.map",
        library = directory("data/centrifuge_db/library")
    threads: 1
    resources:
        mem_mb = 200
    params:
        domain = "archaea,bacteria"
    singularity:
        config["container"]
    log:
        "logs/download_centrifuge_library.log"
    shell:
        """
        centrifuge-download -o {output.library} -m -d {params.domain} refseq > {output.mapping} 2> {log}
        """

rule combine_centrifuge_library_sequences:
    input:
        rules.download_centrifuge_library.output.library
    output:
        temp("data/centrifuge_db/combined.fna")
    threads: 1
    resources:
        mem_mb = 200
    log:
        "logs/combine_centrifuge_library_sequences.log"
    shell:
        """
        cat {input}/*/*.fna > {output} 2> {log}
        """

rule build_centrifuge_db:
    input:
        sequences = "data/centrifuge_db/combined.fna",
        conversion_table = rules.download_centrifuge_library.output.mapping,
        tax_tree = rules.download_centrifuge_taxonomy.output.tax_tree,
        name_table = rules.download_centrifuge_taxonomy.output.name_table
    output:
        expand("data/centrifuge_16s/archaea_bacteria.{db_idx}.cf", db_idx=range(1, 4))
    threads:
        config["build_centrifuge_db"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["build_centrifuge_db"]["memory"]
    params:
        prefix = "data/centrifuge_16s/archaea_bacteria"
    singularity:
        config["container"]
    log:
        "logs/build_centrifuge_db.log"
    shell:
        """
        centrifuge-build -p {threads} --conversion-table {input.conversion_table} \
                 --taxonomy-tree {input.tax_tree} --name-table {input.name_table} \
                 {input.sequences} {params.prefix} &> {log}
        """


rule centrifuge:
    input:
        rules.build_centrifuge_db.output,
        fastq = "data/{run}/filtlong/{sample}.filtered.fastq.gz"
    output:
        report = "data/{run}/centrifuge/centrifuge_report_{run}_{sample}.tsv",
        classification = "data/{run}/centrifuge/centrifuge_classification_{run}_{sample}.tab"
    threads:
        config["centrifuge"]["threads"]
    resources:
        mem_mb = lambda wildcards, attempt: attempt * config["centrifuge"]["memory"]
    params:
        index_prefix = "data/centrifuge_16s/archaea_bacteria"
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
        index_prefix = "data/centrifuge_16s/archaea_bacteria" 
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
          -k 50 \
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
