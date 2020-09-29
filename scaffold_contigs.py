#!//home/eanderson/Virtual_Envs/General3/bin/python3.6

import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_file", help="Input fasta file to scaffold")
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    
    contig_path = args.fasta_file
    scaffolds = set()
    contigs = []

    for record in SeqIO.parse(contig_path, "fasta"):

        record.seq = Seq(str(record.seq), generic_dna)
        record_dict = dict([x.split('=') for x in record.description.split()[1:]])

        for key, val in record_dict.items():
            record_dict[key] = int(val)

        if contigs:
            if record.id in contigs:
                print(f"Contig {record.id} has multiple Fasta entries.", file=sys.stderr)
                sys.exit(1)

            elif int(record.id[1:]) < int(contigs[-1][1:]):
                print(f"Contig {record.id} is smaller than previous contig {contigs[-1]}.", file=sys.stderr)
                print(f"Please input a sorter set of contigs.", file=sys.stderr)
                sys.exit(1)

        contigs.append(record.id)


        if record_dict['SID'] in scaffolds:
            scaff_record.seq = scaff_record.seq + record.seq

        else:
            scaff_record = SeqRecord(record.seq, id = 'SID' + str(record_dict['SID']))
            scaffolds.add(record_dict['SID'])


        if record_dict['GAP'] < 0:
            print(scaff_record.format("fasta"), end="")

        else:
            scaff_record.seq = scaff_record.seq + Seq("N" * record_dict['GAP'], generic_dna)
            
            
if __name__ == "__main__":
    main()


