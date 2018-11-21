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
