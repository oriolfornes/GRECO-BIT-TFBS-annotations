import argparse
import json
import os
import re
import sys

# Append JASPAR-profile-inference to path
jaspar_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                          os.pardir, "JASPAR-profile-inference")
sys.path.append(jaspar_dir)

# Import globals
from __init__ import Jglobals

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an
    {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("tfs", help="TF list (e.g. TAIR, Ensembl, etc.)")
    parser.add_argument("txt", help="UniProt TXT file")
    parser.add_argument("fas", help="UniProt FASTA file")

    taxons = ["fungi", "insects", "nematodes", "plants", "vertebrates"]
    parser.add_argument("tax", choices=taxons, help="Taxon (i.e. %s)" % \
                        ", ".join(taxons), metavar="tax")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Get JASPAR 2 UniProt
    json_file = os.path.join(jaspar_dir, "files", "%s.uniprot.json" % args.tax)
    with open(json_file) as f:
        jaspar2uniprot = json.load(f)

    # Get TFs
    tfs = set()
    for tf in Jglobals.parse_file(args.tfs):
        tfs.add(tf)

    # For each header, sequence...
    fas = {}
    for seq_record in Jglobals.parse_fasta_file(args.fas):
        m = re.search("^\S+\|(\S+)\|(\S+) .+ GN=(.+) PE=\d+ ",
                      seq_record.description)
        if m:
            fas.setdefault(m.group(2), [m.group(1), m.group(3),
                                        str(seq_record.seq)])

    # For each line...
    txt = {}
    for line in Jglobals.parse_file(args.txt):
        if line.startswith("//"):
            if unientry in fas:
                if genes:
                    for gene in genes:
                        txt.setdefault(gene, [])
                        txt[gene].append((uniaccs, unientry, pfams))
                txt.setdefault(fas[unientry][1], [])
                txt[fas[unientry][1]].append((uniaccs, unientry, pfams))
        if line.startswith("ID"):
            m = re.search("^ID\s+(\S+)", line)
            unientry = m.group(1)
            uniaccs = []
            genes = set()
            pfams = set()
        if line.startswith("AC"):
            uniaccs += re.findall("(\S+);", line)
        if line.startswith("GN"):
            m = re.search("OrderedLocusNames=(\S+);", line)
            if m: genes.add(m.group(1))
            m = re.search("ORFNames=(\S+);", line)
            if m: genes.add(m.group(1))
        if line.startswith("DR"):
            m = re.search("^DR\s+Ensembl; \S+; \S+; (\S+).", line)
            if m: genes.add(m.group(1))
            m = re.search("^DR\s+FlyBase; (\S+); \S+.", line)
            if m: genes.add(m.group(1))
            # m = re.search("^DR\s+SGD; \S+; (\S+).", line)
            # if m: genes.add(m.group(1))
            m = re.search("^DR\s+TAIR; locus:\d+; (\S+).", line)
            if m: genes.add(m.group(1))
            m = re.search("^DR\s+WormBase; \S+; \S+; (\S+); \S+.", line)
            if m: genes.add(m.group(1))
            m = re.search("^DR\s+Pfam; PF\d+; (\S+); \S+.", line)
            if m: pfams.add(m.group(1))

    print("Gene Name\tUniProt Accession\tUniProt Entry\tPfam ID\tSequence\tJASPAR")
    # For each tf...
    for tf in sorted(tfs):
        if tf in txt:
            # Initialize
            uniaccs = []
            unientries = []
            families = set()
            sequences = []
            jaspar_motifs = set()
            for uaccs, uentry, pfams in txt[tf]:
                for uniacc in uaccs:
                    if uniacc in jaspar2uniprot:
                        for jaspar_motif in jaspar2uniprot[uniacc][0]:
                            jaspar_motifs.add(jaspar_motif)
                for pfam in pfams:
                    families.add(pfam)
                if fas[uentry][0] not in uniaccs:
                    uniaccs.append(fas[uentry][0])
                    unientries.append(uentry)
                    sequences.append(fas[uentry][2])
                    families.update(pfams)
            print("%s\t%s\t%s\t%s\t%s\t%s" % \
                    (fas[uentry][1], ";".join(uniaccs), ";".join(unientries),
                     ";".join(sorted(families)), ";".join(sequences),
                     ";".join(sorted(jaspar_motifs)))
            )

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()