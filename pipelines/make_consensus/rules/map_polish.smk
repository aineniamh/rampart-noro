rule minimap2_racon0:
    input:
        reads=config["output_path"] + "/binned_{sample}/{analysis_stem}.fastq",
        ref=config["output_path"] + "/binned_{sample}/{analysis_stem}.fasta"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/mapped.paf"
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon1:
    input:
        reads=config["output_path"]+"/binned_{sample}/{analysis_stem}.fastq",
        fasta=config["output_path"] + "/binned_{sample}/{analysis_stem}.fasta",
        paf= config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/mapped.paf"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon1.fasta"
    shell:
        "../racon/build/bin/racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}"

rule minimap2_racon1:
    input:
        reads=config["output_path"]+"/binned_{sample}/{analysis_stem}.fastq",
        ref= config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon1.fasta"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/mapped.racon1.paf"
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon2:
    input:
        reads=config["output_path"]+"/binned_{sample}/{analysis_stem}.fastq",
        fasta= config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon1.fasta",
        paf= config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/mapped.racon1.paf"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon2.fasta"
    shell:
        "../racon/build/bin/racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}"


rule minimap2_racon2:
    input:
        reads=config["output_path"]+"/binned_{sample}/{analysis_stem}.fastq",
        ref= config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon2.fasta"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/mapped.racon2.paf"
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon3:
    input:
        reads=config["output_path"]+"/binned_{sample}/{analysis_stem}.fastq",
        fasta= config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon2.fasta",
        paf= config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/mapped.racon2.paf"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon3.fasta"
    shell:
        "../racon/build/bin/racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}"


rule minimap2_racon3:
    input:
        reads=config["output_path"]+"/binned_{sample}/{analysis_stem}.fastq",
        ref= config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon3.fasta"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/mapped.racon3.paf"
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule racon4:
    input:
        reads=config["output_path"]+"/binned_{sample}/{analysis_stem}.fastq",
        fasta= config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon3.fasta",
        paf= config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/mapped.racon3.paf"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon4.fasta"
    shell:
        "../racon/build/bin/racon --no-trimming -t 1 {input.reads} {input.paf} {input.fasta} > {output}"

rule minimap2_racon4:
    input:
        reads=config["output_path"]+"/binned_{sample}/{analysis_stem}.fastq",
        ref= config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon4.fasta"
    output:
        config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/mapped.racon4.paf"
    shell:
        "minimap2 -x map-ont {input.ref} {input.reads} > {output}"

rule medaka:
    input:
        # basecalls=config["output_path"] + "/binned_{sample}/{analysis_stem}.fastq",
        # draft=config["output_path"] + "/binned_{sample}/{analysis_stem}.fasta"
        basecalls=config["output_path"]+"/binned_{sample}/{analysis_stem}.fastq",
        draft= config["output_path"] + "/binned_{sample}/polishing/{analysis_stem}/racon4.fasta"
    params:
        outdir=config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}"
    output:
        consensus= config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.fasta"
    threads:
        2
    shell:
        "medaka_consensus -i {input.basecalls} -d {input.draft} -o {params.outdir} -t 2"

rule minimap2_medaka:
    input:
        reads=config["output_path"]+"/binned_{sample}/{analysis_stem}.fastq",
        ref= config["output_path"] + "/binned_{sample}/{amplicon}/medaka/consensus.fasta"
    output:
        config["output_path"] + "/binned_{sample}/{amplicon}/medaka/consensus.mapped.sam"
    shell:
        "minimap2 -ax map-ont {input.ref} {input.reads} > {output}"

rule gather_sequences:
    input:
        expand(config["output_path"] + "/binned_{{sample}}/medaka/{analysis_stem}/consensus.fasta", analysis_stem=config["analysis_stem"])
    output:
        config["output_path"]+"/consensus_sequences/{sample}.fasta"
    shell:
        "cat {input} > {output}"