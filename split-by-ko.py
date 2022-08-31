#! /usr/bin/env python
import os
import sys
import argparse
import screed

import pandas as pd
from collections import defaultdict

import sourmash
from sourmash.logging import notify
#from sourmash.lca.lca_utils import make_lineage

def make_outdir(output_dirname):
    if not os.path.exists(output_dirname):
        try:
            os.makedirs(output_dirname)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

def write_fastas(fastaD, infoD, output_dir, fasta_type):
    outdir = os.path.join(output_dir, fasta_type)
    make_outdir(outdir)

    total_ko = len(fastaD.keys())
    print(f'writing {fasta_type} fastas')
    for n, (ident, records) in enumerate(fastaD.items()):
        if n % 10000 == 0:
            progress = float(n)/total_ko * 100
            print(f"writing file for {n}'th ko, {ident} ({progress}% of {fasta_type} ko's).\n")
        # make filename
        this_filename = f"{outdir}/{ident}.fna"
        if fasta_type == "protein":
            this_filename = f"{outdir}/{ident}.faa"
        with open(this_filename, 'w') as out_fa:
            for rc in records:
                out_fa.write(f">{rc.name}\n{rc.sequence}\n")
        # store fromfile info
        infoD[ident][fasta_type] = this_filename
        #info_out.write(f"{this_ident},,{this_filename}\n")

    print(f"{str(n)} {fasta_type} ko's written to fasta files\n")
    fastaD = None # does this help cleanup memory or is it useless?
    return infoD


def split_fasta_dict(fasta):
    koD = defaultdict(list)
    for n, record in enumerate(screed.open(fasta)):
        # get ko
        ko = record.name.rsplit('|ko:', 1)[1]
        # add record to fasta dict
        koD[ko].append(record)
    print(f"parsed {n} records in {fasta}")
    return koD


def main(args):
    notify(f"Splitting fastas into groups by 'ko'. Writing files to '{args.output_dir}' \n")
    output_dir = args.output_dir
    # split fasta by ko:
    infoD = defaultdict(lambda: defaultdict(dict))
    proteinD = split_fasta_dict(args.protein_fasta)
    infoD = write_fastas(proteinD, infoD, output_dir, "protein")
    if args.nucl_fasta:
        nuclD = split_fasta_dict(args.nucl_fasta)
        infoD = write_fastas(proteinD, infoD, output_dir, "nucl")

    # make sourmash sketch fromfile csv
    with open(args.output_fromfile_csv, 'w') as info_out:
        print("now writing 'fromfile' csv")
        info_out.write("name,genome_filename,protein_filename\n") # header needs these three, even if we only have proteins
        for ident, fasta_files in infoD.items():
            prot = fasta_files["protein"]
            nucl = fasta_files.get("nucl", "")
            info_out.write(f"{ident},{nucl},{prot}\n")


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--protein-fasta", default= "kegg_genes_KO.faa")
    p.add_argument("--nucl-fasta", default= "kegg_genes_KO.fna")
    p.add_argument("--output-dir", default="ko_split_fastas")
    p.add_argument("--output-fromfile-csv", default= "ko.fromfile.csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
