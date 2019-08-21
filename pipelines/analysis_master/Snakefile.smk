##### Configuration #####

configfile: "pipelines/analysis_master/config.yaml"

barcodes = config["barcodes"].split(',')

barcode_string = ""
for i in barcodes:
    barcode_string+=" {}".format(i)

sample_string = ""
if config["sample"]:
    sample_string = config["sample"]
else:
    sample_string = "_".join(barcodes)
    config["sample"]=sample_string
    


##### Subworkflows #####
rule all:
    input:
        expand(config["output_path"] + "/consensus_sequences/{sample}.fasta",  sample=config["sample"])

rule bin_to_fastq:
    input:
    params:
        barcodes = config["barcodes"],
        sample= config["sample"],
        path_to_reads= config["input_path"],
        path_to_csv= config["annotated_path"],
        output_path= config["output_path"],
        min_length = config["min_length"],
        max_length = config["max_length"]
    output:
        fastq=config["output_path"] + "/binned_{sample}.fastq",
        csv=config["output_path"] + "/binned_{sample}.csv"
    shell:
        "snakemake --nolock --snakefile pipelines/bin_to_fastq/Snakefile "
        "--config "
        "input_path={params.path_to_reads} "
        "output_path={params.output_path} "
        "annotated_path={params.path_to_csv} "
        "min_length={params.min_length} "
        "max_length={params.max_length} "
        "barcodes={params.barcodes} "
        "sample={params.sample} "

rule assess_sample:
    input:
        fastq=config["output_path"] + "/binned_{sample}.fastq",
        csv= config["output_path"] + "/binned_{sample}.csv",
        refs = config["references_file"],
        config = config["config"]
    params:
        sample = "{sample}",
        output_path= config["output_path"],
        min_reads=config["min_reads"],
        min_pcent=config["min_pcent"]
    output:
        config["output_path"] + "/binned_{sample}/config.yaml"
    shell:
        "snakemake --nolock --snakefile pipelines/assess_sample/Snakefile "
        "--configfile {input.config} "
        "--config "
        "reads={input.fastq} "
        "csv={input.csv} "
        "output_path={params.output_path} "
        "references={input.refs} "
        "config_in={input.config} "
        "config_out={output} "
        "min_reads={params.min_reads} "
        "min_pcent={params.min_pcent} "
        "sample={params.sample} "

rule make_consensus:
    input:
        config=config["output_path"] + "/binned_{sample}/config.yaml"
    params:
        sample = "{sample}",
        output_path= config["output_path"]
    output:
        config["output_path"] + "/consensus_sequences/{sample}.fasta"
    shell:
        "snakemake --nolock --snakefile pipelines/make_consensus/Snakefile "
        "--configfile {input.config} "
        "--config "
        "output_path={params.output_path} "
        "sample={params.sample}"
