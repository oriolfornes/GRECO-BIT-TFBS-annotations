import argparse
# import coreapi
# import json
import os
import pickle
import re
import shutil
import subprocess
import sys

# Append JASPAR-profile-inference to path
jaspar_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                          os.pardir, os.pardir, "JASPAR-profile-inference")
sys.path.insert(0, os.path.realpath(jaspar_dir))

# Import from JASPAR-profile-inference
from __init__ import Jglobals

#-------------#
# Class       #
#-------------#

class TF(object):

    def __init__(self, gene_name):

        self.gene_name = gene_name
        self.uniacc = None
        self.unientry = None
        self.status = None
        self.sequence = None
        self.pfam_id = "Unknown"
        self.cluster_num = None
        # self.orthodb = set()
        self.jaspar_id = []
        self.hocomoco_id = set()
        
        # In vivo
        self.chip_atlas = set()
        self.cistromedb = set()
        self.gtrd = set()
        self.dap_seq = set()
        self.remap = set()

        # In vitro
        self.ht_selex = set()
        self.cisbp = set()
        self.uniprobe = set()
        self.smile_seq = set()

        # Hidden variables (for internal use only)
        self._uniaccs = set()
        self._unientries = set()
        self._sequences = set()
        self._pfam_ids = set()
        self._species = set()

    def __str__(self):

        string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.gene_name,
            self.uniacc,
            ";".join(sorted([i for i in self._uniaccs if i != self.uniacc])),
            self.unientry,
            ";".join(sorted([i for i in self._unientries if i != self.unientry])),
            self.status,
            self.sequence,
            ";".join(sorted([i for i in self._sequences if i != self.sequence])),
            sorted(self._species)[0],
            self.pfam_id,
            self.cluster_num,
            ";".join(sorted([i for i in self._pfam_ids if i != self.pfam_id])),
            # ";".join(sorted(self.orthodb)),
            self.jaspar_id,
            ";".join(sorted(self.hocomoco_id)),
            ";".join(sorted(self.chip_atlas)),
            ";".join(sorted(self.cistromedb)),
            ";".join(sorted(self.gtrd)),
            ";".join(sorted(self.remap)),
            ";".join(sorted(self.dap_seq)),
            ";".join(sorted(self.ht_selex)),
            ";".join(sorted(self.cisbp)),
            ";".join(sorted(self.uniprobe)),
            ";".join(sorted(self.smile_seq))
        )

        return(string)

    def __repr__(self):

        string = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
            self.gene_name,
            self.uniacc,
            ";".join(sorted([i for i in self._uniaccs if i != self.uniacc])),
            self.unientry,
            ";".join(sorted([i for i in self._unientries if i != self.unientry])),
            self.status,
            self.sequence,
            ";".join(sorted([i for i in self._sequences if i != self.sequence])),
            sorted(self._species)[0],
            self.pfam_id,
            self.cluster_num,
            ";".join(sorted([i for i in self._pfam_ids if i != self.pfam_id])),
            # ";".join(sorted(self.orthodb)),
            self.jaspar_id,
            ";".join(sorted(self.hocomoco_id)),
            ";".join(sorted(self.chip_atlas)),
            ";".join(sorted(self.cistromedb)),
            ";".join(sorted(self.gtrd)),
            ";".join(sorted(self.remap)),
            ";".join(sorted(self.dap_seq)),
            ";".join(sorted(self.ht_selex)),
            ";".join(sorted(self.cisbp)),
            ";".join(sorted(self.uniprobe)),
            ";".join(sorted(self.smile_seq))
        )

        return(string)

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an
    {argparse} object.
    """

    parser = argparse.ArgumentParser(description="creates CSV file for GRECO.")

    parser.add_argument("--dummy-dir", default="/tmp/", metavar="DIR",
                        help="dummy directory (default = /tmp/)")
    parser.add_argument("--output-dir", default="./", metavar="DIR",
                        help="output directory (default = ./)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="verbose mode")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Create output dir
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    #-------------#
    # Parse Data  #
    #-------------#

    # Initialize
    tfs = set()
    species = []
    pfam_ids = {}
    hocomoco_file = "./Data/Databases/TFs/HOCOMOCO/HOCOMOCOv11_full_jaspar_format.txt"
    chip_atlas_file = "./Data/Experiments/ChIP-seq.ChIP-Atlas.tsv"
    cistromedb_file = "./Data/Experiments/ChIP-seq.CistromeDB.tsv"
    gtrd_file = "./Data/Experiments/ChIP-seq.GTRD.tsv"
    remap_file = "./Data/Experiments/ChIP-seq.ReMap2020.tsv"
    dap_seq_file = "./Data/Experiments/DAP-seq.PMID:27203113.tsv"
    ht_selex_files = ["./Data/Experiments/HT-SELEX.PMID:23332764.tsv",
                      "./Data/Experiments/HT-SELEX.PMID:28473536.tsv"]
    cisbp_file = "./Data/Experiments/PBM.CIS-BP.tsv"
    uniprobe_file = "./Data/Experiments/PBM.UniPROBE.tsv"
    smile_seq = "./Data/Experiments/SMiLE-seq.PMID:28092692.tsv"
    species_files = ["./Data/Species/Arabidopsis_thaliana_TFs.tab.gz",
                     "./Data/Species/Caenorhabditis_elegans_TFs.tab.gz",
                     "./Data/Species/Homo_sapiens_TFs.tab.gz",
                     "./Data/Species/Drosophila_melanogaster_TFs.tab.gz",
                     "./Data/Species/Mus_musculus_TFs.tab.gz",
                     "./Data/Species/Saccharomyces_cerevisiae_TFs.tab.gz"]

    # # Skip if pickle file already exists
    # pickle_file = os.path.join(args.output_dir, "1.pkl")
    # if not os.path.exists(pickle_file):

    # For each file...
    for file_name in species_files:

        # Get species
        m = re.search("(\w+)_(\w+)_TFs.tab.gz", file_name)
        species.append("%s %s" % (m.group(1), m.group(2)))
        species.append("%s_%s" % (m.group(1), m.group(2)))

        # For each line...
        for line in Jglobals.parse_tsv_file(file_name):

            if line[0] == "Gene Name":
                continue

            # Initialize
            tf = TF(line[0])

            # Get UniProt Accession
            uniaccs = line[1].split(";")
            tf.uniacc = uniaccs[-1]
            tf._uniaccs.update(set(uniaccs))

            # Get UniProt Entry
            unientries = line[2].split(";")
            tf.unientry = unientries[-1]
            tf._unientries.update(set(unientries))

            # Get sequence
            sequences = line[5].split(";")
            tf.sequence = sequences[-1]
            tf._sequences.update(set(sequences))

            # Get Pfam IDs
            try:
                pfam_families = line[4].split(";")
                for pfam_id in pfam_families:
                    tf._pfam_ids.add(pfam_id)
                    pfam_ids.setdefault(pfam_id, set())
                    pfam_ids[pfam_id].add((tf.uniacc, tf.sequence))
            except:
                pass

            # Get species
            tf._species.add(species[-1])
            tf._species.add(species[-2])

            # Get JASPAR ids
            tf.jaspar_id = line[6]

            # # Get orthoDB cluster
            # codec = coreapi.codecs.CoreJSONCodec()
            # for uniacc in tf._uniaccs:
            #     json_file = os.path.join(args.orthodb, "%s.json" % uniacc)
            #     if not os.path.exists(json_file):
            #         client = coreapi.Client()
            #         response = client.get(
            #             "https://www.orthodb.org/search?query=%s&level=2759&species=2759" % uniacc)
            #         json_obj = json.loads(codec.encode(response))
            #         with open(json_file, "w") as j:
            #             j.write(json.dumps(json_obj, sort_keys=True, indent=4, separators=(",", ": ")))
            #     with open(json_file, "r") as j:  
            #         json_obj = json.load(j)
            #         for orthodb in json_obj["data"]:
            #             tf.orthodb.add(orthodb)

            # Add TF to TFs
            tfs.add(tf)

    #     # Write pickle file
    #     handle = Jglobals._get_file_handle(pickle_file, "wb")
    #     pickle.dump(tfs, handle)

    # else:

    #     # Load pickle file
    #     handle = Jglobals._get_file_handle(pickle_file, "rb")
    #     tfs = pickle.load(handle)

    #-------------#
    # HOCOMOCO    #
    #-------------#

    #  678 Homo sapiens
    #  451 Mus musculus

    # SMCA5_MOUSE # Not a TF: SWI/SNF-related matrix-associated actin-dependent regulator of chromatin subfamily A member 5
    # FUBP1_MOUSE # Not a TF: Far upstream element-binding protein 1
    # BRCA1_MOUSE # Not a TF: Breast cancer type 1 susceptibility protein homolog
    # EVI1_MOUSE  # Not a TF: Histone-lysine N-methyltransferase MECOM
    # TAF1_MOUSE  # Not a TF: Transcription initiation factor TFIID subunit 1
    # BRAC_MOUSE  # Not a valid UniProt Entry
    # HLTF_MOUSE  # Not a TF: Helicase-like transcription factor
    # TAF1_HUMAN  # Not a TF: Transcription initiation factor TFIID subunit 1
    # HLTF_HUMAN  # Not a TF: Helicase-like transcription factor
    # BRCA1_HUMAN # Not a TF: Breast cancer type 1 susceptibility protein
    # SMCA5_HUMAN # Not a TF: SWI/SNF-related matrix-associated actin-dependent regulator of chromatin subfamily A member 5
    # ZF64A_HUMAN # Not a valid UniProt Entry
    # BRAC_HUMAN  # Not a valid UniProt Entry
    # EVI1_HUMAN  # Not a TF: Histone-lysine N-methyltransferase MECOM
    # FUBP1_HUMAN # Not a TF: Far upstream element-binding protein 1
    # BPTF_HUMAN  # Not a TF: Nucleosome-remodeling factor subunit BPTF
    # CENPB_HUMAN # Not a TF: Major centromere autoantigen B
    # ZBT48_HUMAN # Not a valid UniProt Entry; should be TZAP_HUMAN

    # Skip if pickle file already exists
    pickle_file = os.path.join(args.output_dir, "2.pkl")
    if not os.path.exists(pickle_file):

        # For each line...
        for line in Jglobals.parse_file(hocomoco_file):

            if line.startswith(">"):

                # Get unientry
                m = re.search("(\w+_(HUMAN|MOUSE))", line)
                unientry = m.group(1)

                # For each TF...
                for tf in sorted(tfs, key=lambda x: x.gene_name):
                    if unientry in tf._unientries:
                        tf.hocomoco_id.add(line[1:])

        # Write pickle file
        handle = Jglobals._get_file_handle(pickle_file, "wb")
        pickle.dump(tfs, handle)

    else:

        # Load pickle file
        handle = Jglobals._get_file_handle(pickle_file, "rb")
        tfs = pickle.load(handle)

    if args.verbose:
        genes = 0
        feats = 0
        print(">HOCOMOCO...")
        for tf in sorted(tfs, key=lambda x: x.gene_name):
            subtotal_feats = len(tf.hocomoco_id)
            if subtotal_feats > 0:
                genes += 1
                feats += subtotal_feats
                print(tf.gene_name, tf.hocomoco_id)
        print("Total genes: %s" % genes)
        print("Total feats: %s" % feats)
        print("//")

    #-------------#
    # GTRD        #
    #-------------#

    #   71 Arabidopsis thaliana
    #  213 Caenorhabditis elegans
    #   11 Danio rerio
    #  249 Drosophila melanogaster
    # 1236 Homo sapiens
    #  513 Mus musculus
    #   12 Rattus norvegicus
    #  137 Saccharomyces cerevisiae
    #   32 Schizosaccharomyces pombe

    # Skip if pickle file already exists
    pickle_file = os.path.join(args.output_dir, "3.pkl")
    if not os.path.exists(pickle_file):

        # For each line...
        for line in Jglobals.parse_tsv_file(gtrd_file):

            # Inialize
            experiment_id = line[0]
            uniacc = line[1]

            # For each TF...
            for tf in sorted(tfs, key=lambda x: x.gene_name):
                if uniacc in tf._uniaccs:
                    tf.gtrd.add(experiment_id)

        # Write pickle file
        handle = Jglobals._get_file_handle(pickle_file, "wb")
        pickle.dump(tfs, handle)

    else:

        # Load pickle file
        handle = Jglobals._get_file_handle(pickle_file, "rb")
        tfs = pickle.load(handle)

    if args.verbose:
        genes = 0
        feats = 0
        print(">GTRD...")
        for tf in sorted(tfs, key=lambda x: x.gene_name):
            subtotal_feats = len(tf.gtrd)
            if subtotal_feats > 0:
                genes += 1
                feats += subtotal_feats
                print(tf.gene_name, tf.gtrd)
        print("Total genes: %s" % genes)
        print("Total feats: %s" % feats)
        print("//")

    #-------------#
    # DAP-seq     #
    #-------------#

    # Skip if pickle file already exists
    pickle_file = os.path.join(args.output_dir, "after_dapseq.pkl")
    if not os.path.exists(pickle_file):

        # # Get DAP-seq SRA runs
        # uniprot = {}
        # # For each file...
        # for file_name in os.listdir(args.uniprot):
        #     if os.path.isfile(os.path.join(args.uniprot, file_name)):
        #         if args.dap_seq.endswith(file_name):
        #             with open(os.path.join(args.uniprot, file_name), "r") as f:
        #                 # For each row...
        #                 for row in csv.reader(f, delimiter="\t"):
        #                     if row[0] == "Entry": continue
        #                     uniprot.setdefault(row[1], row[0])

        # For each line...
        for line in Jglobals.parse_tsv_file(dap_seq_file):

            if line[0] == "AvgSpotLen":
                continue

            species = line[8]
            sra_run = line[9]
            gene = line[15]

            if line[8] not in species:
                continue

            # For each TF...
            for tf in sorted(tfs, key=lambda x: x.gene_name):
                if gene.upper() in tf.gene_name.upper():
                    tf.dap_seq.add(sra_run)

        # Write pickle file
        handle = Jglobals._get_file_handle(pickle_file, "wb")
        pickle.dump(tfs, handle)

    else:

        # Load pickle file
        handle = Jglobals._get_file_handle(pickle_file, "rb")
        tfs = pickle.load(handle)

    if args.verbose:
        genes = 0
        feats = 0
        print(">DAP-seq...")
        for tf in sorted(tfs, key=lambda x: x.gene_name):
            subtotal_feats = len(tf.dap_seq)
            if subtotal_feats > 0:
                genes += 1
                feats += subtotal_feats
                print(tf.gene_name, tf.dap_seq)
        print("Total genes: %s" % genes)
        print("Total feats: %s" % feats)
        print("//")
    exit(0)

    #-------------#
    # HT-SELEX    #
    #-------------#

    # Get HT-SELEX SRA runs
    for file_name in args.ht_selex:
        m = re.search("HT-SELEX.PMID:(\d+).tsv", file_name)
        pmid = int(m.group(1))
        with open(file_name, "r") as f:
            # For each row...
            for row in csv.reader(f, delimiter="\t"):
                if row[0] == "Alias": continue
                m = re.search("^([A-Za-z\d]+)_", row[0])
                if m:
                    gene_name = m.group(1)
                    if pmid == 23332764: sra_run = row[14]
                    if pmid == 28473536: sra_run = row[19]
                    # For each TF...
                    for tf in sorted(tfs, key=lambda x: x.gene_name):
                        # HT-SELEX data only available for human and mouse
                        if "Homo sapiens" not in tf._species and "Mus musculus" not in tf._species: continue
                        if tf.gene_name == gene_name:
                            tf.ht_selex.add(sra_run)

#    for tf in sorted(tfs, key=lambda x: x.gene_name):
#        if len(tf.ht_selex) > 0: print(tf.gene_name)
#    exit(0)

    #-------------#
    # CIS-BP      #
    #-------------#

    # Get CIS-BP ids
    species = re.search("(\w+_\w+).tab", args.tfang)
    valid_matrices = set()
    for file_name in os.listdir(args.cisbp_dir):
        if os.path.isfile(os.path.join(args.cisbp_dir, file_name)):
            m = re.search("^(M\d{4}_\d.\d{2})", file_name)
            if m:
                valid_matrices.add(m.group(1))
    with open(args.cisbp, "r") as f:
        # For each row...
        for row in csv.reader(f, delimiter="\t"):
            if row[0] == "TF_ID": continue
            if row[3] in valid_matrices and row[8] == "D":
                matrix_id = row[3]
                gene_name = row[6]
                # For each TF...
                for tf in sorted(tfs, key=lambda x: x.gene_name):
                    if tf.gene_name.upper() == gene_name.upper() and row[7] in tf._species:
                        tf.cisbp.add(matrix_id)

#    for tf in sorted(tfs, key=lambda x: x.gene_name):
#        if len(tf.cisbp) > 0: print(tf.gene_name)
#    exit(0)

    # Get UniPROBE ids
    with open(args.uniprobe, "r") as f:
        # For each row...
        for row in csv.reader(f, delimiter="\t"):
            if row[0] == "Protein": continue
            gene_name = row[0]
            uniprobe_id = row[1]
            # For each TF...
            for tf in sorted(tfs, key=lambda x: x.gene_name):
                if tf.gene_name.upper() == gene_name.upper() and row[2] in tf._species:
                    tf.uniprobe.add(uniprobe_id)

#    for tf in sorted(tfs, key=lambda x: x.gene_name):
#        if len(tf.uniprobe) > 0: print(tf.gene_name)
#    exit(0)

    # Get SMiLE-seq SRA runs
    if args.smile_seq:
        synonyms = {
            "CEBPb": "CEBPB",
            "cFOS": "FOS",
            "cFOSL2": "FOSL2",
            "cJUN": "JUN",
            "PPARa": "PPARA",
            "PPARg": "PPARG",
            "RXRa": "RXRA",
            "RXRg": "RXRG",            
        }
        with open(args.smile_seq, "r") as f:
            # For each row...
            for row in csv.reader(f, delimiter="\t"):
                if row[0] == "Assay_Type": continue
                sra_run = row[13]
                for gene_name in row[15].split("_")[0].split("-"):
                    if gene_name in synonyms:
                        gene_name = synonyms[gene_name]
                    # For each TF...
                    for tf in sorted(tfs, key=lambda x: x.gene_name):
                        if tf.gene_name.upper() == gene_name.upper() and row[12] in tf._species:
                            tf.smile_seq.add(sra_run)


    #-------------#
    # Cluster TFs #
    #-------------#

    # Initialize
    pid = os.getpid()
    fasta_file = os.path.join(args.dummy_dir, "sequences_%s.fa" % pid)
    prefix = os.path.join(args.dummy_dir, "sequences_%s" % pid)
    cluster_all_seqs = "%s_all_seqs.fasta" % prefix
    cluster_clusters = "%s_cluster.tsv" % prefix
    cluster_rep_seqs = "%s_rep_seq.fasta" % prefix

    # For each Pfam ID...
    for pfam_id in sorted(pfam_ids, key=lambda x: len(pfam_ids[x]), reverse=True):

        # Initialize
        cluster_num = 0
        clusters = {}

        # For each TF...
        for uniacc, sequence in pfam_ids[pfam_id]:
            Jglobals.write(fasta_file, ">%s\n%s" % (uniacc, sequence))

        # Default parameters for kClust (PMID:23945046):
        # Three similarity criteria are used to decide if a query is added to a cluster:
        #
        #   (1) The sequence similarity score from the 4-mer-based dynamic programming algorithm
        #       is larger than a minimum BLOSUM62 score per column (default 1.12 half bits, which
        #       corresponds to a sequence identity of 30%, see Additional file 1: Figure S2);
        #   (2) the alignment achieves an E-value less than a defined threshold (default value 1E-3); and
        #   (3) the alignment covers at least 80% of the residues of the representative sequence.
        #
        # This criterion ensures that clusters contain sequences with nearly identical domain composition.
        #
        # However, with "-c 0.5", clusters match the TFCLASS classification better.
        opts = "--min-seq-id 0.3 -c 0.5"
        cmd = "mmseqs easy-cluster %s %s %s %s" % (fasta_file, prefix,
                                                   prefix, opts)
        process = subprocess.run(cmd, shell=True, check=True)

        # For each line...
        for line in Jglobals.parse_file(cluster_clusters):

            # Get cluster
            cluster_id, uniacc = line.split("\t")
            if cluster_id not in clusters:
                cluster_num += 1
                clusters.setdefault(cluster_id, cluster_num)

            # For each TF...
            for tf in tfs:
                if uniacc in tf._uniaccs and tf.pfam_id == "Unknown":
                    tf.pfam_id = pfam_id
                tf.cluster_num = clusters[cluster_id]

        # Remove files
        shutil.rmtree(prefix)
        os.remove(cluster_all_seqs)
        os.remove(cluster_clusters)
        os.remove(cluster_rep_seqs)

#    for tf in sorted(tfs, key=lambda x: x.gene_name):
#        if len(tf.smile_seq) > 0: print(tf.gene_name)
#    exit(0)

    # For each TF...
    for tf in sorted(tfs, key=lambda x: x.gene_name):
        print(tf)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()