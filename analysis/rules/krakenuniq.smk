rule download_krakenuniq_db:
    output:
        "data/krakenuniq_db/taxonomy/names.dmp",
        "data/krakenuniq_db/taxonomy/nodes.dmp"
    threads:
        cluster_config["download_krakenuniq_db"]["nCPUs"]
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
