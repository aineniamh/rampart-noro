# assess sample

A snakemake pipeline that takes in the filtered fastq containing all reads for a particular barcode/sample and assesses how many parallel analyses should be run, binning the reads by reference if a significant portion of the sample hits that reference. Follows on from the bin_to_fastq snakemake pipeline.


```
snakemake --snakefile pipelines/assess_sample/Snakefile \
--configfile pipelines/assess_sample/config.yaml \
--config \
input_path=../artic-polio/data/test_rampart_binned  \
output_path=../artic-polio/data/test_rampart_binned \
min_pcent=2 min_reads=100 sample=sample_19_20 \
config=pipelines/assess_sample/config.yaml \
references=references.fasta
```
