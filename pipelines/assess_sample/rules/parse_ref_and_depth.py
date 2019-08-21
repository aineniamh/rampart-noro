import argparse
from Bio import SeqIO
import collections
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

parser = argparse.ArgumentParser(description='Parse mappings, add to headings and create report.')

parser.add_argument("--csv", action="store", type=str, dest="report")
parser.add_argument("--reads", action="store", type=str, dest="reads")
parser.add_argument("--references", action="store", type=str, dest="references")

parser.add_argument("--sample", action="store", type=str, dest="sample")
parser.add_argument("--min_reads", action="store", type=int, dest="min_reads")
parser.add_argument("--min_pcent", action="store", type=float, dest="min_pcent")

parser.add_argument("--output_path", action="store", type=str, dest="output_path")

parser.add_argument("--config_in", action="store", type=str, dest="config_in")
parser.add_argument("--config_out", action="store", type=str, dest="config_out")


args = parser.parse_args()

def make_ref_dict(references):
    refs = {}
    for record in SeqIO.parse(references,"fasta"):
        refs[record.id]=record.seq
    return refs

ref_dict = make_ref_dict(str(args.references))

unknown=False
unmapped_count=0
coord_unmapped = 0
record_count=0

report = pd.read_csv(args.report)
ref_count  = report['best_reference'].value_counts()
fig, ax = plt.subplots(figsize=(15,11))
sns.barplot(ref_count.index, ref_count.values, alpha=0.8)
plt.title('Reference profile of sample {}'.format(args.sample))
plt.ylabel('Read Count', fontsize=12)
plt.xlabel('Reference', fontsize=12)
plt.xticks(rotation=20)
fig.savefig(args.output_path + "/reference_count.pdf")


total = len(report)
refs = []
for i,x in zip(list(ref_count.index), list(ref_count.values)):
    pcent = 100*(x/total)
    if x>args.min_reads and pcent > args.min_pcent:
        refs.append(i)

with open(args.config_out,"w") as new_config: #file to write genotype information
    with open(args.config_in, "r") as f:
        for l in f:
            l=l.rstrip('\n')
            # if l.startswith("barcodes"):
            #     new_config.write("barcode: {}\n".format(args.sample))
            if not l.startswith("config"):
                new_config.write(l + '\n')
    new_config.write("analysis_stem:\n")
    for ref in refs:
        new_config.write("  - {}\n".format(ref))

for ref in refs:
    with open(args.output_path + "/" + ref+ ".fasta","w") as fw:
        fw.write(">{}\n{}\n".format(ref, ref_dict[ref]))

    filtered_df = report.loc[(report["best_reference"]==ref)]
    read_names = list(filtered_df["read_name"].values)
    new_file = args.output_path + "/" + ref + ".fastq"
    with open(new_file,"w") as fw:
        records = []
        for record in SeqIO.parse(args.reads,"fastq"):
            if record.id in read_names:
                records.append(record)

        SeqIO.write(records, fw, "fastq")
