rule minimap2_medaka:
    input:
        reads=output_dir+"/binned/barcode_{barcode}.fastq",
        ref= output_dir + "/medaka/{barcode}/consensus.fasta"
    output:
        output_dir + "/medaka/{barcode}/consensus.mapped.bam"
    shell:
        "minimap2 -ax map-ont {input.ref} {input.reads} | samtools view -b - > {output}"

#rule variants

#rule phasing

#rule minor