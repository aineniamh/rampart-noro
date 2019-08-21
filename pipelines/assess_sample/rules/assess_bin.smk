rule assess_sample:
    input:
        reads= config["output_path"] + "/binned_{sample}.fastq",
        csv= config["output_path"] + "/binned_{sample}.csv",
        refs = config["references"],
        config = config["config"]
    params:
        sample = "{sample}",
        output_path = config["output_path"] + "/binned_{sample}",
        min_reads = config["min_reads"],
        min_pcent = config["min_pcent"],
        path_to_script = workflow.current_basedir
    output:
        fig = config["output_path"] + "/binned_{sample}/reference_count.pdf",
        new_config =  config["output_path"] + "/binned_{sample}/config.yaml"
    shell:
        "python {params.path_to_script}/parse_ref_and_depth.py "
        "--reads {input.reads} "
        "--csv {input.csv} "
        "--output_path {params.output_path} "
        "--references {input.refs} "
        "--config_in {input.config} "
        "--config_out {output.new_config} "
        "--min_reads {params.min_reads} "
        "--min_pcent {params.min_pcent} "
        "--sample {params.sample} "
        






