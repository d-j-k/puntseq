build_bracken_db:
    input:
        rules.build_kraken2_db.output
    output:
        "data/kraken2_db/database35mers.kmer_distrib"
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
