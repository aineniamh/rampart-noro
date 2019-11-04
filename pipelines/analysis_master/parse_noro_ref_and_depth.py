import argparse
from Bio import SeqIO
import collections
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

parser = argparse.ArgumentParser(description='Parse mappings, add to headings and create report.')

parser.add_argument("--csv", action="store", type=str, dest="report")
parser.add_argument("--reads", action="store", type=str, dest="reads")
parser.add_argument("--references", action="store", type=str, dest="references")
parser.add_argument("--amplicons", action="store", type=str, dest="amplicons")
parser.add_argument("--reference_options", action="store", type=str, dest="reference_options")


parser.add_argument("--sample", action="store", type=str, dest="sample")
parser.add_argument("--min_reads", action="store", type=int, dest="min_reads")
parser.add_argument("--min_pcent", action="store", type=float, dest="min_pcent")

parser.add_argument("--output_path", action="store", type=str, dest="output_path")


args = parser.parse_args()

def parse_reference_file(references):
    #returns a dict of dicts containing reference header information
    #key is seq id i.e. first field of the header string
    ref_info = defaultdict(dict)
    for record in SeqIO.parse(references,"fasta"):
        tokens = record.description.split(' ')
        for i in tokens:
            try:
                info = i.split('=')
                ref_info[record.id][info[0]]=info[1]
                # ref_info['*'][info[0]]='NA'
            except:
                pass
    return ref_info

def parse_reference_options(reference_options):
    #returns a dict of {key:list} pairs containing  
    # csv_header_to_be : list of corresponding read_header_group with optional coordinates
    # e.g. "genogroup" : [["genogroup"]]
    # e.g. "loc_genotype" : [["POL_genogroup",0,5000],["VP_genogroup",5000,7000]]
    columns = reference_options.split(";")
    ref_options = defaultdict(list)
    for i in columns:
        k,v = i.rstrip(']').split("[")
        v = [i.split(':') for i in v.split(",")]
        new_v =[]
        for sublist in v:
            if len(sublist)==3:
                new_v.append([sublist[0],int(sublist[1]),int(sublist[2])])
            else:
                new_v.append(sublist)
        ref_options[k]=new_v
        
    return ref_options, ','+','.join(ref_options.keys())


def make_amp_dict(amplicons):
    amp_dict = collections.defaultdict(list)
    with open(amplicons, "r") as f:
        for l in f:
            l=l.rstrip("\n")
            tokens= l.split(",")
            accession,reference,amplicon= tokens[0].split('|')
            coords=tokens[1].split(':')
            if amplicon in ["Amp1","Amp2","Amp3"]:
                amp_dict["Amp123"].append(int(coords[0])).append(int(coords[1]))
            if amplicon in ["Amp4"]:
                amp_dict["Amp4"].append(int(coords[0])).append(int(coords[1]))
            if amplicon in ["Amp5"]:
                amp_dict["Amp5"].append(int(coords[0])).append(int(coords[1]))
    amp_coords = {}
    for amp in amp_dict:
        s= sorted(amp_dict[amp])
        start,end = s[0],s[-1]
        amp_coords[amp] = (start,end)
    return amp_coords

def check_overlap(coords1,coords2):
    list1 = list(range(coords1[0],coords1[1]))
    list2 = list(range(coords2[0],coords2[1]))
    overlap = set(list1).intersection(list2)
    if overlap:
        return True, len(overlap)
    else:
        return False, 0 

amp_dict = make_amp_dict(str(args.amplicons))
ref_dict = make_ref_dict(str(args.references))

def parse_line(line):

    values = {}

    tokens = line.rstrip('\n').split(',')
    values["read_name"], values["read_len"] = tokens[:2]
    values["ref_hit"], values["ref_len"], values["coord_start"], values["coord_end"], values["matches"], values["aln_block_len"] = tokens[5:11]
    try:
        values["ref_options"] = tokens[11:]
    except:
        pass

    return values

def parse_csv(report, reads, reference_options, reference_info):
    #This function parses the input paf file 
    #and outputs a csv report containing information relevant for RAMPART and barcode information
    # read_name,read_len,start_time,barcode,best_reference,start_coords,end_coords,ref_len,matches,aln_block_len,ref_option1,ref_option2
    counts = {
        "unmapped": 0,
        "total": 0,
        "Amp123": collections.Counter(),
        "Amp4": collections.Counter(),
        "Amp5": collections.Counter()
    }

    with open(str(report),"r") as f:
        for line in f:
            counts["total"]+=1
            mapping = parse_line(line)
            if mapping["ref_hit"] in ['*','?','none']:
                counts["unmapped"]+=1
            else:
                overlap_list = []
                for i in amp_dict:
                    overlap, length = check_overlap((int(mapping["coord_start"]),int(mapping["coord_end"])),amp_dict[i])
                    if overlap:
                        overlap_list.append(i, length)
                amplicon = sorted(overlap_list, key = lambda x : x[1], reverse=True)[0][0]
                counts[amplicon][mapping["ref_hit"]]+=1
    



unknown=False
unmapped_count=0
coord_unmapped = 0
record_count=0

report = pd.read_csv(args.report)
report["ref_stem"]=report["best_reference"].str.split("_").str[0]
detailed_ref_count = report['best_reference'].value_counts()

ref_count  = report['ref_stem'].value_counts()
if len(ref_count) > 1:
    fig, ax = plt.subplots(figsize=(15,11))
    sns.barplot(ref_count.index, ref_count.values, alpha=0.8)
    plt.title('Reference profile of sample {}'.format(args.sample))
    plt.ylabel('Read Count', fontsize=12)
    plt.xlabel('Reference', fontsize=12)
    plt.xticks(rotation=20)
    fig.savefig(args.output_path + "/reference_count.pdf")

detail_dict= {}
for i,x in zip(list(detailed_ref_count.index), list(detailed_ref_count.values)):
    stem = i.split("_")[0]
    if stem not in detail_dict:
        detail_dict[stem] = i


total = len(report)
refs = []
for i,x in zip(list(ref_count.index), list(ref_count.values)):
    pcent = 100*(x/total)
    if x>args.min_reads and pcent > args.min_pcent:
        if i != '*':

            refs.append(i.split('_')[0])

print(",".join(refs))

for ref in refs:
    with open(args.output_path + "/" + ref+ ".fasta","w") as fw:
        
        fw.write(">{} detail={}\n{}\n".format(ref, detail_dict[ref], ref_dict[detail_dict[ref]]))

    filtered_df = report.loc[(report["ref_stem"]==ref)]
    read_names = list(filtered_df["read_name"].values)
    new_file = args.output_path + "/" + ref + ".fastq"
    with open(new_file,"w") as fw:
        records = []
        for record in SeqIO.parse(args.reads,"fastq"):
            if record.id in read_names:
                records.append(record)

        SeqIO.write(records, fw, "fastq")

none_df = report.loc[(report["best_reference"]=='*')]
read_names = list(none_df["read_name"].values)
new_file = args.output_path + "/no_hit.fastq"
with open(new_file,"w") as fw:
    records = []
    for record in SeqIO.parse(args.reads,"fastq"):
        if record.id in read_names:
            records.append(record)

    SeqIO.write(records, fw, "fastq")