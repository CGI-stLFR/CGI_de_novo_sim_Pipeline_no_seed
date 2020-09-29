#!/home/eanderson/Virtual_Envs/General3/bin/python3

'''
Calculate the phase-block N50 from phased contigs
'''

import argparse
import sys

import numpy as np
from Bio import SeqIO

def process_args(args):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile", help="The input FASTA file")
    parser.add_argument("--outfile", default=sys.stdout, type=argparse.FileType('w'), help="The output file")
    parser.add_argument("--n_metrics", nargs="+", type=float, default=[0.5], help="NXX phase-block metrics to calculate")
    return parser.parse_args(args)

def get_tuple(contig_info, idx, step):
    while idx >= 0 and idx < len(contig_info):
        contig_tuple = contig_info[idx]
        if contig_tuple[1] != '-1':
            return contig_tuple
        idx += step
    return ('', '', 0)


def main(args):
    scaffold_dict = {}
    contig_info = []
    # Find the length of all heterozygous phase-blocks
    for record in SeqIO.parse(open(args.infile), "fasta"):
        record_dict = dict([x.split('=') for x in record.description.split('\t')[1:]])
        scaffold, pid = record_dict["SID"], record_dict["PID"]
        contig_len = len(record)
        contig_info.append((scaffold, pid, contig_len))
        # Handle the PID=-1 contigs in the second pass
        if pid == "-1":
            continue
        # Add the contig length to it's respective PID
        pid_dict = scaffold_dict.setdefault(scaffold, {})
        if pid not in pid_dict:
            pid_dict[pid] = 0
        pid_dict[pid] += len(record)
    
    # Add the homozygous contigs to the larger adjecent phase-block
    for i, contig_tuple in enumerate(contig_info):
        if contig_tuple[1] != "-1":
            continue
        previous_tuple = get_tuple(contig_info, i - 1, -1)
        next_tuple = get_tuple(contig_info, i + 1, 1)

        found_pid = None
        if next_tuple[0] != contig_tuple[0] and previous_tuple[0] != contig_tuple[0]:
            print(f"Contig: {contig_tuple}; has no adjecent contigs. Next: {next_tuple}. Previous: {previous_tuple}", file=sys.stderr)
            continue
        elif next_tuple[0] != contig_tuple[0]:
            found_pid = previous_tuple[1]
        elif previous_tuple[0] != contig_tuple[0]:
            found_pid = next_tuple[1]
        else: # All have the same SID - match the larger phase_block
            next_size = scaffold_dict[contig_tuple[0]][next_tuple[1]]
            previous_size = scaffold_dict[contig_tuple[0]][previous_tuple[1]]
            if next_size > previous_size:
                found_pid = next_tuple[1]
            else:
                found_pid = previous_tuple[1]
        scaffold_dict[contig_tuple[0]][found_pid] += contig_tuple[2]
    
    # Calculate the phase-block NXX
    phase_blocks = np.array(list(reversed(sorted([x for phase_dict in scaffold_dict.values() for x in phase_dict.values()]))))

    print(f"Minimum length: {phase_blocks.min()}", file=args.outfile)
    print(f"Maximum length: {phase_blocks.max()}", file=args.outfile)
    print(f"Mean length: {phase_blocks.mean()}", file=args.outfile)
    print(f"Median length: {np.median(phase_blocks)}", file=args.outfile)
    print(f"Number: {str(len(phase_blocks))}", file=args.outfile)
    total_bases = np.sum(phase_blocks)
    print(f"Total bases: {total_bases}", file=args.outfile)
    for cutoff in args.n_metrics:
        cutoff_bases = total_bases * cutoff
        current_bases = 0
        for length in phase_blocks:
            current_bases += length
            if current_bases > cutoff_bases:
                break
        print(f"N{int(cutoff * 100)} phase-block length: {length}")

if __name__ == "__main__":
    args = process_args(None)
    main(args)

