rule map_to_cns:
    input:
        cns = rules.medaka.output,
        reads = rules.files.params.reads
    output:
        paf = config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.mapped.paf"
    shell:
        "minimap2 -x map-ont {input.cns} {input.reads} > {output}"

rule mask_low_coverage_regions:
    input:
        cns = rules.medaka.output,
        paf = config["output_path"] + "/binned_{sample}/medaka/{analysis_stem}/consensus.mapped.paf"
    params:
        path_to_script = workflow.current_basedir
    output:
        config["output_path"] +"/binned_{sample}/consensus_sequences/{analysis_stem}.fasta"
    shell:
        """
        python {params.path_to_script}/mask_low_coverage.py \
        --cns {input.cns} \
        --paf {input.paf} \
        --min_coverage 100 \
        --masked_cns {output}
        """

rule cat_stems:
    input:
        expand(config["output_path"] +"/binned_{{sample}}/consensus_sequences/{analysis_stem}.fasta",analysis_stem=config["analysis_stem"])
    output:
        config["output_path"] + "/consensus_sequences/{sample}.fasta"
    shell:
        "cat {input} > {output}"
