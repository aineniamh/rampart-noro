from Bio import SeqIO
from Bio import AlignIO
import sys
import argparse

parser = argparse.ArgumentParser(description='Clean up coding consensus.')
parser.add_argument("--consensus", action="store", type=str, dest="consensus")
parser.add_argument("--alignment_with_ref", action="store", type=str, dest="alignment_with_ref")
parser.add_argument("--output_seq", action="store", type=str, dest="output_seq")
parser.add_argument("--polish_round", action="store", type=str, dest="round")
args = parser.parse_args()

round_name = ''
if args.round:
    round_name = f" round_name={args.round}"


def find_gaps(aln):
    gap_dict = {}
    alignment = AlignIO.read(aln, "fasta")

    print(f"Reading in {alignment}.\nLooking for gaps in coding sequence.")
    for i in range(len(alignment[0])):

        col = alignment[:, i]
        if len(set(col)) >1:
            print(alignment[:, i])
            if '-' in col:
                print(f"Gap at position {i}: {col}")
                if col[0]=='-':
                    gap_dict[i] = ''
                else:
                    gap_dict[i] = col.rstrip('-')+'N'
    return(gap_dict)

#the rule is to replace a gap in the query with 'N' and to force delete a base that causes a gap in the reference
with open(args.output_seq, "w") as fw:

    gap_dict = find_gaps(args.alignment_with_ref)
    
    for record in SeqIO.parse(args.consensus, "fasta"):
        

        new_seq = list(record.seq)

        for key in gap_dict:
            try:
                new_seq[key]= gap_dict[key]
            except IndexError:
                new_seq.append(gap_dict[key])
        new_seq = ''.join(new_seq).upper()
        

        fw.write(f">{record.id}{round_name} length={len(new_seq)}\n{new_seq}\n")

