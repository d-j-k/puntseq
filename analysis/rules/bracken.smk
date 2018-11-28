build_bracken_db:
    input:
        rules.build_kraken2_db.output
    threads:
        cluster_config["build_bracken_db"]["nCPUs"]
    resources:
        mem_mb = cluster_config["build_bracken_db"]["memory"]
    params:
        kraken2_db = "data/kraken2_db"
    singularity:
        config["container"]
    log:
        "logs/build_bracken_db.log"
    shell:
        """
        bracken-build -d {params.kraken2_db} \
          -threads {threads} \
          -l ${READ_LEN} 2> {log}
        """
