#!/usr/bin/env python

import argparse
# import coreapi
import json
import os
import re
import sys

# Get root path
root_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        os.pardir, os.pardir)

# Append JASPAR-profile-inference to path
jaspar_dir = os.path.join(root_dir, "JASPAR-profile-inference")
sys.path.insert(1, os.path.realpath(jaspar_dir))

# Import from JASPAR-profile-inference
from __init__ import Jglobals

#-------------#
# Class       #
#-------------#

class CISBP(object):

    def __init__(self, id, name, species, geneid, family, evidence):

        self.id = id
        self.name = name
        self.species = species
        self.geneid = geneid
        self._family = family
        self.evidence = evidence

    @property
    def family(self):

        if self._family == "AP-2":
            return("AP2")

        return(self._family)

    def __str__(self):

        string = "{}\t{}\t{}\t{}\t{}\t{}".format(
            self.id,
            self.name,
            self.species,
            self.geneid,
            self.family,
            self.evidence
        )

        return(string)

    def __repr__(self):

        string = "{}\t{}\t{}\t{}\t{}\t{}".format(
            self.id,
            self.name,
            self.species,
            self.geneid,
            self.family,
            self.evidence
        )

        return(string)

# #-------------#
# # Functions   #
# #-------------#

# def uniacc_to_OrthoDB_cluster(uniacc):

#     # Initialize
#     codec = coreapi.codecs.CoreJSONCodec()
#     client = coreapi.Client()
#     orthodb_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
#                                ".OrthoDB")

#     if not os.path.isdir(orthodb_dir):
#         os.mkdir(orthodb_dir)

#     try:

#         # Skip if already JSON data already fetched
#         json_file = os.path.join(orthodb_dir, ".%s.json" % uniacc)
#         if not os.path.exists(json_file):

#             # Get response
#             url = "https://www.orthodb.org/search?query=%s&level=2759&species=2759" % \
#                 uniacc
#             response = client.get(url)
#             json_obj = json.loads(codec.encode(response))
#             with open(json_file, "w") as j:
#                 j.write(json.dumps(json_obj))

#         with open(json_file, "r") as j:  
#             json_obj = json.load(j)
#             for orthodb_cluster in json_obj["data"]:
#                 yield(orthodb_cluster)

#     except:
#         pass

#-------------#
# Main        #
#-------------#

# Initialize
species = {
    "Arabidopsis_thaliana": "plants",
    "Caenorhabditis_elegans": "nematodes",
    "Drosophila_melanogaster": "insects",
    "Homo_sapiens": "vertebrates",
    "Mus_musculus": "vertebrates",
    "Saccharomyces_cerevisiae": "fungi"
}

for s in sorted(species):

    # Initialize
    lines = []
    ss = s.replace("_", " ")

    # Get JASPAR 2 UniProt
    json_file = os.path.join(jaspar_dir, "files", "%s.uniprot.json" % \
                             species[s])
    with open(json_file) as f:
        jaspar2uniprot = json.load(f)

    # Get TFs
    tfs = set()
    tfs_file = os.path.join(root_dir, "Data", "Databases", "CisBP-2.0",
                            "%s.csv" % s)
    for line in Jglobals.parse_csv_file(tfs_file):
        if line[0] == "ID":
            continue
        cisbp = CISBP(line[0], line[1], line[2], line[3], line[4], line[5])
        tfs.add(cisbp)

    # For each SeqRecord...
    fas = {}
    fas_file = os.path.join(root_dir, "Data", "Databases", "UniProt",
                            "%s_TFs.fa.gz" % s)
    for seq_record in Jglobals.parse_fasta_file(fas_file):
        m = re.search("^\S+\|(\S+)\|(\S+) .+ OX=\d+ GN=(.+) PE=\d+ ",
                        seq_record.description)
        if m:
            fas.setdefault(m.group(2), [m.group(1), m.group(3),
                                        str(seq_record.seq)])

    # For each line...
    txt = {}
    txt_file = os.path.join(root_dir, "Data", "Databases", "UniProt",
                            "%s_TFs.txt.gz" % s)
    for line in Jglobals.parse_file(txt_file):
        if line.startswith("//"):
            if unientry in fas:
                for gene in genes:
                    txt.setdefault(gene, {})
                    txt[gene].setdefault(reviewed, [])
                    # txt[gene][reviewed].append((uniaccs, unientry, orthodb))
                    txt[gene][reviewed].append((uniaccs, unientry, geneid))
        if line.startswith("ID"):
            m = re.search("^ID\s+(\S+)\s+(Reviewed|Unreviewed);", line)
            unientry = m.group(1)
            reviewed = m.group(2)
            uniaccs = []
            geneid = ""
            genes = set()
            # orthodb = set()
        if line.startswith("AC"):
            uniaccs += re.findall("(\S+);", line)
        if line.startswith("GN"):
            if species[s] == "fungi":
                m = re.search("OrderedLocusNames=(\S+);", line)
                if m:
                    genes.add(m.group(1))
        if line.startswith("DR"):
            m = re.search("^DR\s+GeneID; (\d+); -.", line)
            if m:
                geneid = m.group(1)
                continue
            if species[s] == "insects":
                m = re.search("^DR\s+FlyBase; (\S+); \S+.", line)
                if m:
                    genes.add(m.group(1))
            elif species[s] == "nematodes": 
                m = re.search("^DR\s+WormBase; \S+; \S+; (\S+); \S+.", line)
                if m:
                    genes.add(m.group(1))
            elif species[s] == "plants":
                m = re.search("^DR\s+Araport; (\S+);", line)
                if m:
                    genes.add(m.group(1))
            elif species[s] == "vertebrates":
                m = re.search("^DR\s+Ensembl; ENS\w+; ENS\w+; (ENS\w+)", line)
                if m:
                    genes.add(m.group(1))
            # m = re.search("^DR\s+OrthoDB; (\d+at2759)")
            # if m:
            #     orthodb.add(m.group(1))

    # Initialize
    lines = set()

    # For each CIS-BP object...
    for cisbp in sorted(tfs, key=lambda x: x.geneid):

        if cisbp.geneid not in txt:
            continue

        # Initialize
        uniaccs = []
        unientries = []
        geneids = set()
        family = cisbp.family
        sequences = []
        jaspar_motifs = set()
        for reviewed in ["Reviewed", "Unreviewed"]:
            if reviewed not in txt[cisbp.geneid]:
                continue
            for uaccs, uentry , geneid in txt[cisbp.geneid][reviewed]:
                geneids.add(geneid)
                for uniacc in uaccs:
                    if uniacc in jaspar2uniprot:
                        for jaspar_motif in jaspar2uniprot[uniacc][0]:
                            jaspar_motifs.add(jaspar_motif)
                if fas[uentry][0] not in uniaccs:
                    uniaccs.append(fas[uentry][0])
                    unientries.append(uentry)
                    sequences.append(fas[uentry][2])

            break
        lines.add("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (\
                  fas[uentry][1], ss, ";".join(uniaccs), ";".join(unientries),
                  ";".join(sorted(map(str, geneids))), ";".join(sequences),
                  reviewed, family, ";".join(sorted(jaspar_motifs))))

    # Fix species-specific cases
    if s == "Arabidopsis_thaliana":
        # Q9MAB7_ARATH
        gene = "T4P13.29"
        uniacc = "Q9MAB7"
        unientry = "Q9MAB7_ARATH"
        seq = ["MEEQQEFEKPIFEFRPKKLRTSRVSRNLMKKTGFRESMNHYEELSCNYGLRENPKKTQKS",
               "LLDHRFCTRRRRNKKILIRCKECGKGFLYEKCLFNHLQVTHSEESTRRSLFCRFSIVQRR",
               "KRSKRVSRYKKILPRFSVSSSSCTMFPVSVDDGGLLEVAESLILLSMSGGKFVNGLEHFG",
               "KALGSTQRKFEDGLLRNEQRLVGEALDSNPEKKLVGISRASVGTSKELSGYLANKKGRED",
               "DELGQQKQAGARILREETDNEQKLVRQETAFEDSVSGFEMNIEHRCGLCHKVFSTYQTLG",
               "GHQTFHRMRNKSKSQTKRCREESIEAEAGINNGSVTLTISAEFAEGCLGQYML"]
        reviewed = "Unreviewed"
        family = "C2H2 ZF"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT5G14340
        gene = "MYB40"
        uniacc = "F4K6R6"
        unientry = "F4K6R6_ARATH"
        seq = ["MCDLLACLRKISMGRKPCCDKIGLKRGPWTIEEDHRLMNFILNNGIHCWRIVPKLAGLLR",
               "CGKSCRLRWINYLRPDLKRGGFTDAEEDRIMELHSQLGNRWSKIASHFSGRTDNEIKNHW",
               "NTKIKKKMKHLGLDPATHKPMNDITHQTDPNQDKKPNMCSTINEGEEIKDQTPKDDVITE",
               "TTKTLMLSDNDEELVAKNCKILCAEEVDLESLFETQCNEISSSSFSSLCSNISRSESSSY",
               "LAEDSISLEQWDLDMTDPFVPWDLFANLDDNLFLL"]
        reviewed = "Unreviewed"
        family = "Myb/SANT"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT4G04555
        gene = "At4g04555"
        uniacc = "A0A1P8B772"
        unientry = "A0A1P8B772_ARATH"
        seq = ["MRGRIVRSYTRSKVPHWRWTDDLNLLFIQVVELLGGERRATPKVILDFMDVKNLPISHVK",
               "SHLQMYRNKKKEESRKERRMMREMSRRQSQQYIQIYERYNWILVRR"]
        reviewed = "Unreviewed"
        family = "Myb/SANT"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT4G04605
        gene = "At4g04605"
        uniacc = "A0A1P8B8R6"
        unientry = "A0A1P8B8R6_ARATH"
        seq = ["MKNPIVRSYVRSKVPRLRWNSDLHNSFVQAVEQLGGENRATPKMVLQLMDVRGLTISHVK",
               "SHLQMYRGMKLEESMQEEILAKRSVRVTGQVTWWQFQHYLHNYQRLQGNANFFQNQQRQE",
               "AGIYEKIITFGQSSNETKEVPSDYFTCKSLYETSRAPQNKTRDDDENGEAVVAIEDDGNV",
               "VDEHDDEDANDDEVTLSLNSIKVSKSEEELSLELTLGLKA"]
        reviewed = "Unreviewed"
        family = "Myb/SANT"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AtMYB84
        gene = "AtMYB84"
        uniacc = "O49746"
        unientry = "O49746_ARATH"
        seq = ["MGRAPCCDKANVKKGPWSPEEDAKLKSYIENSGTGGNWIALPQKIGLKRCGKSCRLRWLN",
               "YLRPNIKHGGFSEEEENIICSLYLTIGSRWSIIAAQLPGRTDNDIKNYWNTRLKKKLINK",
               "QRKELQEACMEQQEMMVMMKRQHQQQQIQTSFMMRQDQTMFTWPLHHHNVQVPALFRIKP",
               "TRFATKKMLSQCSSRTWSRSKIKNWRKQTSSSSRFNDNAFDHLSFSQLLLDPNHNHLGSG",
               "EGFSMNSILSANTNSPLLNTSNDNQWFGNFQAETVNLFSGASTSTSADQSTISWEDISSL",
               "VYSDSKQFF"]
        reviewed = "Unreviewed"
        family = "Myb/SANT"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT1G63490
        gene = "At1g63490"
        uniacc = "Q94BQ7"
        unientry = "Q94BQ7_ARATH"
        seq = ["MYGNDLDTSVYGSGFPRIGDQRPESVEADIWDEYCGSPWNLNNMPKLKGSMLQAIRHNIN",
               "GVTVPWLYLGMLFSSFCWHFEDHCFYSVNYLHWGEAKCWYGIPGSAASAFEKVMRKTLPD",
               "LFDAQPDLLFQLVTMLSPTVLQENKVPVYTVLQEPGNFVITFPKSFHAGFNFGLNCAEAV",
               "NFATADWLPYGGSGAELYRLYRKPSVISHEELLCVVAKGNCCNNEGSIHLKKELLRIYSK",
               "EKTWREQLWKSGILRSSPMFVPECADSVGIEEDPTCIICQQFLHLSAIVCNCRPSVFACL",
               "EHWKHLCECEPTKLRLEYRYTLAELDMMVQEVEKFGGCKTQETKISQRPSSGTKRSIALN",
               "KKQEGMQVSQARPADKWLLRASKVLDAAFSSVEYATLLKESEQFLWAGSEMDRVRDVTKS",
               "LNKAKIWAEAVSDCLSKVEGEVNDDSMKVHLEFIDELLRVNPVPCFNSGYLKLKDYAEEA",
               "RKLSEKIDSALSSSPTITQLELLHSEVSRSPISLKKHEILSKKISSAKMLAKRAKRYLTD",
               "AKPPGIEMDALFKLNSEMLELHVQLPETEGILDLVKKSESARDKSNKVLTGSLSLENVEE",
               "LLHEFDSFSINVPELNILRQYHVDTLSWISRFNDVMVDVREGKDQRKLISDLSSLLRDGA",
               "SLGIQVEGLPLVEVELKKASCREKARTVYTARKSLDFIEQLLSEAVILHIEEEEIFVEIS",
               "GILSTARCWEERASTILENETQMYELKDLVRMSVNIDAVLPTLQGIENTISSAETWLQKS",
               "EPFLSATSSMASSPCSMLELPVLKDLVTQAKLLNVQLQEPRILETLLLNCERWQCDNHQL",
               "LQETEDLLDNAKIDDGTHSNILPKIMDLITRVDSARRSGLALGLNFDELPKLRTASLKLG",
               "WCCKTITLSSSSPTSELLEDVGKPSLQHIQQHLKEGQTLEILPEEYYLGKRLMELKDTGL",
               "EWAKRARKVVTDSGALALEDVFELISEGENLPVHAEQELQSLRARSMLHCICLKPYNSRS",
               "MVSCSQCGEWYHTYCLKLHWRPKAYVCSACCPLAETTPQIDPARATEPERPSLNQRRTRM",
               "VATDAAVNDLKWKTRKHIKRTTKRSPQVHILPWFFT"]
        reviewed = "Unreviewed"
        family = "ARID/BRIGHT"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT1G78635
        gene = "At1g78635"
        uniacc = "A0A178WKG1"
        unientry = "A0A178WKG1_ARATH"
        seq = ["MPLFVLSKSSAKKSQTFMSQEQQGTSPLNTDEALKSHHSFSQNPNYIESNLSRRYQTNFS",
               "LGSQRISGQFLMIPNPGFVDTSSMTNTHDTSLANFQTVPTTTNPKIILKEEKENEEEQQY",
               "WTIKKELTKSDACQRLTLSKSSVEEHILKHLLPEDSQKIDKGKPGITVKVYDNDTDTEHE",
               "LCLAFQCSYVLKNGWVKTFVKRRGLEEGDMIGLFWECSTSKLHFSVLFRVNAKAPAAEEK",
               "EAY"]
        reviewed = "Unreviewed"
        family = "B3"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT2G13985
        gene = "At2g13985"
        uniacc = "A0A1P8AY78"
        unientry = "A0A1P8AY78_ARATH"
        seq = ["MCGEMILIEEKGRSWTVDLKRKNSCPTTYIKRGWRSFCNANGFRAGSFFTFKLFQRGGTL",
               "GLRLFHRELEEEAMPIECLSTEPDRNQNERSSQIWKDSSSPSQNRFVTIL"]
        reviewed = "Unreviewed"
        family = "B3"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT4G31685
        gene = "At4g31685"
        uniacc = "A0A1P8B795"
        unientry = "A0A1P8B795_ARATH"
        seq = ["MMLKVDSRGAVYIIGGNDWKSFCAANEVGAGESLALELIQGGVSPLLKFCSKMEQPSFKA",
               "EDGRHKRARVQNRSQETDKGAETSRASTMGPKLEITEKGEPSRASTMRPKVEIREKIAET",
               "GEPSRASNKSSGIEGNLQHTKPCSVKTDQLAKVKESVVDTLTSIGRFQAELETMKRKLED",
               "SLQELN"]
        reviewed = "Unreviewed"
        family = "B3"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT5G01305
        gene = "BHLH140"
        uniacc = "Q9M041"
        unientry = "BH140_ARATH"
        seq = ["MDDFNLRSENPNSSSTTSSSSSSFHRHKSETGNTKRSRSTSTLSTDPQSVAARDRRHRIS",
               "DRFKILQSMVPGGAKMDTVSMLDEAISYVKFLKAQIWYHQNMLLFINDHETTSSCTYSPG",
               "AGEFGPKLFGYDDDYAPIMDTYSQGVPLTVADSKYTPWFGSVDDEQEHVTYFKYRRATRH",
               "ALRGHCNCIIGETEEFADQREKMEVQIEESGKNQTSPESIEADKAKQIVVLLIGPPGSGK",
               "STFCDTAMRSSHRPWSRICQDIVNNGKAGTKAQCLKMATDSLREGKSVFIDRCNLDREQR",
               "SEFIKLGGPEFEVHAVVLELPAQVCISRSVKRTGHEGNLQGGRAAAVVNKMLQSKELPKV",
               "NEGFSRIMFCYSDADVDNAVNMYNKLGPMDTLPSGCFGEKKLDTKSQPGIMKFFKKVSAL",
               "PASSSNEATNTTRKADEMTANVRVSPVKLGSADIVPTLAFPSISTADFQFDLEKASDIIV",
               "EKAEEFLSKLGTARLVLVDLSRGSKILSLVKAKASQKNIDSAKFFTFVGDITKLRSEGGL",
               "HCNVIANATNWRLKPGGGGVNAAIFKAAGPDLETATRVRANTLLPGKAVVVPLPSTCPLH",
               "NAEGITHVIHVLGPNMNPNRPDNLNNDYTKGCKTLREAYTSLFEGFLSVVQDQSKLPKRS",
               "SQTAVSDSGEDIKEDSERNKKYKGSQDKAVTNNLESESLEDTRGSGKKMSKGWNTWALAL",
               "HSIAMHPERHENVVLEYLDNIVVINDQYPKARKHVLVLARQESLDGLEDVRKENLQLLQE",
               "MHNVGLKWVDRFQNEDASLIFRLGYHSVPSMRQLHLHVISQDFNSDSLKNKKHWNSFTTS",
               "FFRDSVDVLEEVNSQGKANVASEDLLKGELRCNRCRSAHPNIPKLKSHVRSCHSQFPDHL",
               "LQNNRLVARAET"]
        reviewed = "Unreviewed"
        family = "bHLH"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT1G35490
        gene = "At1g35490"
        uniacc = "C0SUZ3"
        unientry = "C0SUZ3_ARATH"
        seq = ["MENLRRLSNPGNFGFIGRSQSRVPQKNMENNISPPNNMHHHSASLDDLFTEDQPAWLDEL",
               "LSEPASPKINKGHRRSASDTAAYLNSALMPSKENHVAGSSWQFQNYDLWQSNSYEQHNKL",
               "GWDFSTANGTNIQRNMSCGALNMSSKPIEKHVSKMKEGTSTKPDGPRSKTDSKRIKHQNA",
               "HRARLRRLEYISDLERTIQVLQVEGCEMSSAIHYLDQQLLMLSMENRALKQRMDSLAEIQ",
               "KLKHVEQQLLEREIGNLQFRRHQQQPQQNQKQVQAIQNRYNKYQPPVTQEPDAQFAALAI"]
        reviewed = "Unreviewed"
        family = "bZIP"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT5G46010
        gene = "At5g46010"
        uniacc = "Q9FNM1"
        unientry = "Q9FNM1_ARATH"
        seq = ["MSSSSNKNSNKDHDDHHNDQGQHNTIHTPSFMHTFEISNISSPSSPVVSDPKPEWKPNQH",
               "QAQILEELFIGGTVNPSLTSIKQITIKLQSYGEEVDDADVYKWFHNRKYSRKPKLVSYFM",
               "FL"]
        reviewed = "Unreviewed"
        family = "Homeodomain"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT4G37435
        gene = "At4g37435"
        uniacc = "A0A1P8B5J4"
        unientry = "A0A1P8B5J4_ARATH"
        seq = ["MGRVKIQIKRINDRQQRNIAFAKRKNGLLKKAYELFVLCNVPVALILFSPSGKLFVFYAK",
               "ERPEEIICKFLARAVRTRLPHEPNIQRLVMQMRSETKCSDESLDDGLRYLIQPKEIEDEI",
               "EEVNARLAEVEKQLERFLEIPDQWESLAELKRREEDFQKTLDIVQSRKRDLKPGNFVSGI",
               "KKFLRI"]
        reviewed = "Unreviewed"
        family = "MADS box"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT5G26865
        gene = "F2P16.19"
        uniacc = "O04632"
        unientry = "O04632_ARATH"
        seq = ["MREDTMFKKALELSTLCDIEVCVILYSRDGELIKTWPEDQSKVRDMAERFSKLHERERRK",
               "KRTNLSLFLRKKILDNSKLSEKVLEMKDSLESGLRVLQDKLLLLQPEKNQTELGQIPVIN",
               "NGQNHW"]
        reviewed = "Unreviewed"
        family = "MADS box"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT3G57980
        gene = "T10K17.190"
        uniacc = "Q9M2Q0"
        unientry = "Q9M2Q0_ARATH"
        seq = ["MEELLLACAVHRHGTDSWDSVASEVHKQNSTFRTLTAIDCRHKYNDLKRRFSRNLVSPGS",
               "ADEETLAAEISSVPWLEELRKLRVDELRREVERYDLSISSLQLKVKTLEDEREKSLKTEN",
               "SDLDRIAETKENHTESGNNSGVPVTELKNSPDPNDNSPGTGSENTNRAVKIAEPVDEEPN",
               "RIGGEDNDEKPAREDSGRGSCESVAKESDRAEPKREGNDSPELVESMDESKGEEDTKETS",
               "DGQSSASFPRKETVDQDQPDNKDQSLTVNKIFVESQPLSDFIEILQSHPIGSHFSRRLET",
               "QETSDYYRIIRQHIDFEMIRSRVEEGYYKTARTKFFRDLLLLINNVRVFYGEPSPEFNAA",
               "KQLYQLIKKQMSFKIPKQTLPPPKEDALVTSKEEVKVSSLKPTLSVPIIACRKRSSLAVR",
               "SPASVTETLKKKTRVVPTVDEKQVSEEEEGRPSDKDEKPIVSKKMARGAAPSTAKKVGSR",
               "NVKTSLNAGISNRGRSPNGSSVLKKSVQQKKGINTSGGSKKQSAASFLKRMKGVSSSETV",
               "VETVKAESSNGKRGAEQRKSNSKSEKVDAVKLPAGQKRLTGKRPTIEKGSPTKKNSGVAS",
               "KRGTASLMAKRDSETSEKETGSSTRPKKRSKR"]
        reviewed = "Unreviewed"
        family = "Myb/SANT"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT5G03780
        gene = "TRFL10"
        uniacc = "C0SVN2"
        unientry = "C0SVN2_ARATH"
        seq = ["MIEKKLSKMKMRTKSIPVRSSRKSLDSLLDDEANENTQCDERDVPHKKRRCLGTSETTDR",
               "GGSVEPLLDLDACIVCEVSDERVSRCCGVDCLLSFHGECLYADLGSTSSSSSSSSEDVSN",
               "PFCPYCWLKIVALKSKTLREKTVEAEKAVCKYLDKEMKSRDEDITLSGDEIGNQEQSTDI",
               "VSDHELQGEKDGCSSKPDADQGKVGTGKVIDEVGASEKVATEKFQDAEDDETAKDQGTRI",
               "LNTGAGKKREVSSFLSMQESFSAKEQDQVQQNEKRRRRGLKIIDSDISSKGSSNERNGED",
               "VTEQVTSSVQVTSPSGRMRNQQATTKVVAKSKTVRDISFFKMDQRRRLLWTYEEEEMLKV",
               "GVEKFAAEANKNMPWRKILEMGEKVFHETRTPADLKDKWRSMVKIMNKNEQGSTLTPTAM"]
        reviewed = "Unreviewed"
        family = "Myb/SANT"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT3G01600
        gene = "NAC044"
        uniacc = "F4J4R5"
        unientry = "F4J4R5_ARATH"
        seq = ["MARAWIVDGRGIAAKVKNASLSSALQIQDCGAHIKCPNCTYRIDNSNVLIPWPGLPKGVK",
               "FEPTDEDIIEFLEAKCGIGGSEPHVLIEEFIRPVTEDVGINYTHPQNLPGANKDGVSVFF",
               "FHKTVQAYGTGQRKRRKITPTLVNDEPVRWHKTGRTKPVMLSGVQRGCKKIMVLYKSARK",
               "GTKPEKSNWVLHQYHLGTEGKEIGDYVVSKITYQQQKLGENPDEGESSSGVRGGPTTPKT",
               "NTPTPPSLVDGVAGDEEAFDDLKMFDPFFEELDSIPEAALGKMWSKKARMDEEFVVNLSE",
               "DNLICDESMEASSLWENQVLPNPSLGTVGDFDGFSISDLENADLGTPPDFLTLASQESLL",
               "NWIGWL"]
        reviewed = "Unreviewed"
        family = "NAC/NAM"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT3G46565
        gene = "At3g46565"
        uniacc = "A0A1I9LRJ5"
        unientry = "A0A1I9LRJ5_ARATH"
        seq = ["MKDVYAKEPWLLDHSNNSFFKEDEWYYFSTRTQISEKKIGRGEITGDNNNGIGRGNWRVN",
               "EKEYIIDKDIGDIIGIKKNLTDKSGN"]
        reviewed = "Unreviewed"
        family = "NAC/NAM"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # AT5G04395
        gene = "At5g04395"
        uniacc = "A0A1P8BGF2"
        unientry = "A0A1P8BGF2_ARATH"
        seq = ["MARDSKKAGGALPAPATAMGNAAAETSLPPGFRFHPSDEELISYYLKKKVQGKPMRYDEI",
               "GEVDICKLEPWDLAVIIFLCFLALDFRYVLKTRDKEWFFFSALDKKTRTGTSMSRATKQG",
               "YWKVTGTDGKIRQGGDGKVTIGTMKTLVFHRGRSPNGLGTDWVMNEYHLAKNDEGVPVRY",
               "MFSTRSF"]
        reviewed = "Unreviewed"
        family = "NAC/NAM"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
    elif s == "Caenorhabditis_elegans":
        # HLH27
        gene = "hlh-27"
        uniacc = "Q18056"
        unientry = "Q18056_CAEEL"
        seq = ["MPKVIPSSMSDYRSVPYNQTPKSGSERKRRNITNELINECKTIVQKSEEEHISQEVVLFR",
               "IVKLVTGVNLESNFSSNDLSESTRRKFDTESERRKVKTEREKIRRKKQDDCYAELKFFIL",
               "NKQMGSYEQRLKLERITILEIIIDYIKHNSDLLYPETIPQILPLLAGKSTATCENKENEK",
               "PKTRMEVKDLFPRLTFQEVQESPTSTSPLLTFPCIPMIPTTQFNVLSNYNTVPSIFSAPL",
               "RFILPSLQILTPETSDEEENEETVDISN"]
        reviewed = "Unreviewed"
        family = "bHLH"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
    elif s == "Drosophila_melanogaster":
        # FBgn0000413
        gene = "da"
        uniacc = "P11420"
        unientry = "DA_DROME"
        seq = ["MATSDDEPMHLYEVFQNCFNKIANKQPTGTVGADRGGGGGYHSPYGSLGVENGMYPSDFN",
               "SMHDTVNGGNNRYANASTVDQYFDSAAAGSGGAWCQPQMSSANSYMGQSAYQNSGPLSGH",
               "SIDQQQQQVHQADGLGMGGGGGGGVGADGMHCPVTTGLPPISSFRPTSGGIGGPGAGQQA",
               "PVNVNVNPPAVFNSPQAHNHNHTVQAQHSALSTAGPLGHHSLNHTPHAHSHTLPLPHALP",
               "HGHTLPHPHHSQQNSPAVQSSDAFSGAGASVKVAGAGNSSAAALRQQMYMPADQSISSFG",
               "SNPSTPVNSPPPLTQSVVGGGGEPSVSGGSGWGHSVLNGGPSSSYASEMVPVSSLHTMAS",
               "VFQGVRMEERLDDALNVLRNHCEPEMLAGVNQSLASIDNIDALTSFVPNSPSHLGSGGNS",
               "GSVSNTSNAALVHEVLALGAAAAAGTSGQSVGGAGSLASLKLDRSASTSLPKQTKKRKEH",
               "TAISNSVPAGVSTTSSLTSLDISDTKPTSSIESSNSGLQQHSQGKGTKRPRRYCSSADED",
               "DDAEPAVKAIREKERRQANNARERIRIRDINEALKELGRMCMTHLKSDKPQTKLGILNMA",
               "VEVIMTLEQQVRERNLNPKAACLKRREEEKAEDGPKLSAQHHMIPQPQQVGGTPGSSYHS",
               "QPAQLVPPSSQTISTMTISLPVNQANNGLPPHLQQQQQQQSQLGHAQLPQ"]
        reviewed = "Reviewed"
        family = "bHLH"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # FBgn0029173
        gene = "prg"
        uniacc = "Q9VLK8"
        unientry = "Q9VLK8_DROME"
        seq = ["MHTSDISMTQDQDVSTCRLCHHNTDPNSLNIFDDTVQFCKDVSIAEVSKSLWSVQYDRNE",
               "CLSELICSRCLEILEEAFELRKGMQEREQSLQEQLKEMIKDHPKHRPGLNGNPGVFVPEE",
               "GCIIVEVDPENLAESSEEEFALGSDGEYENYDDDDEEEEEDYDEEDEEDGQNGEDVDMPL",
               "GMDAAQMAAQQSVANNANTTEARPKRAFLCQYCDLGFTLPAECQEHELAAHDPNAPYCCN",
               "FCNIKLVTRPALISHIKTLHDPDRPYVCAHCRKGFVRRSDLKKHTIVHTGVRPFTCNVCS",
               "KSFSRNTNLTKHMRIHSGVKPFVCQQCPRSFQTAVEMMRHTRSHGEVKAFQCGRCPYSFS",
               "RRDKLIAHQQVHTRRDMEQQQQMGLIPPMEGDLQQQALQAKQKAAAQTKNSRYYHCDVCD",
               "RTFQRERDLQRHQALHMDSLFACKTCNQGFNRREQLQRHELEAHGPSFTCGICCISFLHQ",
               "IELENHLKVHQLQHKMAQRAQEAAILPLKMAEKAPVAMTAPLVQDPQLVRPSAAELSFYS",
               "NMIPTMNLGFYSETRPEE"]
        reviewed = "Unreviewed"
        family = "C2H2 ZF"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # FBgn0261613
        gene = "Oaz"
        uniacc = "A1Z9R4"
        unientry = "ZN423_DROME"
        seq = ["MMALKMLYRGPSSRLENLIEKIQATKEITNNDMYSTHTSSSYSPSISDGTMTPNSHHLIG",
               "APTAAGQEDHPTEGKINGGADGEDLPKPKRLPHFHHHHHHHYHHQQALKIANKLRKINKE",
               "AKMGATAGGGATGAASKFDKLTGEGIKSRGDGSYQCQFCEKTFPRLGYLKHHVQSHAEHL",
               "PFKCEYCSKLFKHKRSRDRHKKLHTNERNYKCPHCEAAFSRSDHLKIHMKTHDIQKPFQC",
               "SMCNRGYNTAAALTSHMQKHKKNAAILAAGGNPNALNYSPRSTGSASASVSSNGSLQKRR",
               "YALALASDSSPSRMDFPKRSRSNHVGGTTTTATPTPLLRCSYCPKVTEFSSLEQLNAHLQ",
               "SVHEQPQTQAVKTPVQEGEGFQLSCEYCTMKFGNIAGLFQHMRSTHMDRLSSPNSYYEHF",
               "NRLATAGTFSPRLALDLPKIKPDLGSPERESRPAEDDLPTDLSNNKRRPLTPNPQAQTPL",
               "APPSAPPGVFFCNQCNAGLPDFESFRNHLKSHIAEGMQLVCPHCGMSLPEQSEFERHVVG",
               "HFLITGSEFNCSSSCGKSFAKSEDLQQHLLSEHVLTLLKCSLCSELCESRMAMQLHLACA",
               "HSQETKLLRCSACLELFRSDAEFHVHVKTRHQLGGHPTLGATSSAPTNPLQCMFCRAVCS",
               "SELEMHFHLAAHARQFRCPSCPETFHVEFLLDRHMQSQHGGVKDKEANSPNMGSLYVNAL",
               "LPPLAAAAAAAAATNNNSSIIDYNVAFKGLFGGASGGAGSGGGGAQSGGAPPSANKFYSP",
               "LQVDTNALKAQTSPHPALMYGLSQRYLMEMYAAKSTSPSGNEGVGNSQPPAPQATAPPPP",
               "PNASTATFSCGMCERQDLRSEAELHSHRKLAHNLKTGVSLRCAYCAGNFKSRAELEQHMK",
               "SCHNSTGKHKCLICDEVFPSPAILAEHKLQHSKVGQSGKCSHCGQPLEDVAAFRAHLSEH",
               "GSDGASLPLACICCRQTLHSEFELSLHAKFHTKSSSSGGSLQEPVCALCLEPLPDATEGP",
               "AKLCDKCCRKHNLNGKRGKHSEPATSLPAPPSAFVENRCNLCKMILPHAQKLQEHLVEHT",
               "FAGTEQRGFNCYICSAVFTAPGGLLNHMGEHGAHSRPYDCNLCPEKFFFRAELEHHQRGH",
               "ELRPQARPPAAKVEVPSIRNTSPGQSPVRSPTIVKQELYETDTVESAGVEDEPENHPDEE",
               "EYIEVEQMPHETRPSGIGSQLERSTSSA"]
        reviewed = "Reviewed"
        family = "C2H2 ZF"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # FBgn0022720
        gene = "zf30C"
        uniacc = "Q9VL91"
        unientry = "Q9VL91_DROME"
        seq = ["METESMAPSSSAAAAATKLLVEITTDGIPKLHCPVCNKALVSLAGYVKHVKKHQPPGGFE",
               "CRHCDARFCHEEELTQHAKDEHGVTGAVAGQERKPFVCEKCGAEYKYQEAYRRHCRTKCG",
               "EEKLPREESRPMECKCCYTRFSSASNLSKHRRSRPDTCGQPEYDSPGSSDGMKKHKAFRK",
               "KDRNRDSDDEDTTSEESEDSDDDIPLASRLKTKLKQESQNSDSGDECPDFEPNNSEDDAD",
               "ASGFQLPPPAMVKVEAFDEEDFEYQDASMYVKTESTDIFSNEKDKLLDVLLNEGDGLKPF",
               "ESLKVEQGAGILDEIAAVPLVEVAEEDVLELRGHQMEKPPGPRKRGRPPKEKIPVVKRKY",
               "RKRNAPRSPSPDDSSGTPKRVAKISKKELKERLKMINKMEKSWKCPHCVKIYHIRKPYEK",
               "HLRDDHKLNEAEMKEIFKDVDVHAKLDEEVFKCPICSKIYLVEKRLVTHMKVHGEDGKLT",
               "FKCPCYCNLFFATKEQATEHARAQHKELLYCEKCDKYMTGHDSLKNHERNFHSKKEPRSQ",
               "QRNLICDKCGKKFTGRTSLSDHVRSDCGRLPLYGCSVCGKHLSTAGILKTHMLLHKADTP",
               "YQCDKCGKTFKVKAQYKSHLKTRHTDYKPYKCHLCPKEYPYRESLLTHMTVHTGIKRFLC",
               "NNCGKRFTCISNLQAHRKVHADTCGQLPLNAKATQYMGVQRGKLLMGAKPEAGMEYEETK",
               "TLIAQDVIDRDMPMAQELNFPSDGSAPLATVPLNYASTHLVPHIVLQTMQAEARRME"]
        reviewed = "Unreviewed"
        family = "C2H2 ZF"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # BRC1_DROME
        gene = "br"
        uniacc = "Q01295"
        unientry = "BRC1_DROME"
        seq = ["MDDTQHFCLRWNNYQSSITSAFENLRDDEAFVDVTLACEGRSIKAHRVVLSACSPYFREL",
               "LKSTPCKHPVILLQDVNFMDLHALVEFIYHGEVNVHQKSLQSFLKTAEVLRVSGLTQQQA",
               "EDTHSHLAQIQNLANSGGRTPLNTHTQSLPHPHHGSLHDDGGSSTLFSRQGAGSPPPTAV",
               "PSLPSHINNQLLKRMAMMHRSSAAAAAEETSHAFKRLRGSDNSLPLSGAVGSGSNNNSPD",
               "LPPLHARSASPQQTPADFSTIKHHNNNNTPPLKEEKRNGPTGNGNSGNGNGNGNGASNGN",
               "GISISDKLGSLTPSPLARAGADDVKSEPMDMVCSNNNANANDEHSNDSTGEHDANRSSSG",
               "DGGKGSLSSGNDEEIGDGLASHHAAPQFIMSPAENKMFHAAAFNFPNIDPSALLGLNTQL",
               "QQSGDLAVSPQGKLTTNATTTTTTINNSITNNNNNNNNNNYDYSLPTKNSNSQKTPSPTT",
               "TTLTTPTTTTPTRPTAITSASGICGLNLSTFAANGSSSGGSNGGLSMTALLPQQQQQQQQ",
               "QHQMSQQQQQQQQQQQQQGNSSSGQQQQPNGILACSTPKANTPTTTQQQMYAAVMAAAAS",
               "ASASTSGSANSSLNNSNSTLNTSGGLNNSASGGDDFRCNPCNKNLSSLTRLKRHIQNVHM",
               "RPTKEPVCNICKRVYSSLNSLRNHKSIYHRNLKQPKQEPGVGATQAAANSFYHQQHQQQQ",
               "LNHHSSS"]
        reviewed = "Reviewed"
        family = "C2H2 ZF"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # FBgn0028993
        gene = "scro"
        uniacc = "Q9NCG1"
        unientry = "Q9NCG1_DROME"
        seq = ["NLSSIHHLQNLHSQHQSTLFNSNHSTPFSVTDILSPIEESYRKLELNGNPPSPFRSNSSS",
               "SSINSPGTLTTSTMANPYAMGTLYHSPGVQTYCGPTDNLSLAGHYTDMRNSASWYGSTAN",
               "DPRFAISRLMSSSASGTMSHMGNMSGLAACSVSDSKPLQFPLAQRRKRRVLFTQAQVYEL",
               "ERRFKQQRYLSAPEREHLASLIHLTPTQVKIWFQNHRYKCKRQAKRKAMAEQNQHNQPAS",
               "SPRRVAVPVLVKDGKPCSGNNSSSQSQQHGTNSTSAGNNTGSANNGNANSGIVSVTANVS",
               "GGLNLITGDAPNSHSPDTSSSLLASYGTVGGSNVAMLQQPCNNTLMSNSLAMAYRNQNNF",
               "IFNGHQQQCGGYLPLQGRAW"]
        reviewed = "Unreviewed"
        family = "Homeodomain"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # FBgn0005630
        gene = "lola"
        uniacc = "P42283"
        unientry = "LOLA1_DROME"
        seq = ["MDDDQQFCLRWNNHQSTLISVFDTLLENETLVDCTLAAEGKFLKAHKVVLSACSPYFATL",
               "LQEQYDKHPIFILKDVKYQELRAMMDYMYRGEVNISQDQLAALLKAAESLQIKGLSDNRT",
               "GGGVAPKPESSGHHRGGKLSGAYTLEQTKRARLATGGAMDTSGDVSGSREGSSSPSRRRR",
               "KVRRRSMENDAHDNSNSSVLQAAASNQSILQQTGAGLAVSALVTTQLSSGPAAGTSSQAS",
               "STQQQQPLTSTNVTKKTESAKLTSSTAAPASGASASAAVQQAHLHQQQAQTTSDAINTEN",
               "VQAQSQGGAQGVQGDDEDIDEGSAVGGPNSATGPNPASASASAVHAGVVVKQLASVVDKS",
               "SSNHKHKIKDNSVSSVGSEMVIEPKAEYDDDAHDENVEDLTLDEEDMTMEELDQTAGTSQ",
               "GGEGSSQTYATWQHDRSQDELGLMAQDAQQRDPQDLSRKENTAPDVASTAEIQRSFQRSI",
               "LNGKQRDEQKIQLPGSRRKRLSVTEVSDMLFEFYKTKSAKVPKAEQPHRQVSPTSGEILD",
               "PSTISAIAVYGTASETASKNLNADEVMRVQNATATRVVGAAAGAAASFHPRPKYTLKTAA",
               "SSTEHTTAIPTSVLVANSAAALTPKPQAAVIAEALMRNGLHNFQQQLRAQEILRQQTPHR",
               "RIKEENDVEIAGGDITPTKILENLLRKQQERDLRHSECENEPGYSTEDDEEGRYHAFDDI",
               "HLMEQSGGKFGNNSGMGMFNANAHGGSASSILDAHQAFRNLEFTLSDYGGSSSNGSTTSP",
               "NGIGLDGEPVYECRHCGKKYRWKSTLRRHENVECGGKEPSHQCPYCPYKSKQRGNLGVHV",
               "RKHHTDLPQLPSKRRSKYSMNRENGMSGSMSDDSQGKLIIDFNGKGELETK"]
        reviewed = "Reviewed"
        family = "C2H2 ZF"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # FBgn0011764
        gene = "Dsp1"
        uniacc = "Q24537"
        unientry = "HMG2_DROME"
        seq = ["MEHFHQIQQTIQHYQQQLAAQQQQQVQQQQLQQHQVVVQQNQQQAHQNSSNTTAGVGTQQ",
               "LFTYKMASSFPNPATTMAQVVATSNAAGTTGYDYRLNMAQAAAAAAVPGSQWWYSAANQG",
               "QVDANTAAQLQHQQQQQQQQQQQQQQQHQQQQQMQQQQQQQNVINSASPMSRVKADAKPR",
               "GRMTAYAYFVQTCREEHKKKHPDETVIFAEFSRKCAERWKTMVDKEKKRFHEMAEKDKQR",
               "YEAEMQNYVPPKGAVVGRGKKRKQIKDPNAPKRSLSAFFWFCNDERNKVKALNPEFGVGD",
               "IAKELGRKWSDVDPEVKQKYESMAERDKARYEREMTEYKTSGKIAMSAPSMQASMQAQAQ",
               "KAALLAAAAQQQHQQLEEQHDDDDGDGDDDENQ"]
        reviewed = "Reviewed"
        family = "Sox"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # FBgn0000099
        gene = "ap"
        uniacc = "P29673"
        unientry = "APTE_DROME"
        seq = ["MGVCTEERPVMHWQQSARFLGPGAREKSPTPPVAHQGSNQCGSAAGANNNHPLFRACSSS",
               "SCPDICDHSTKPFGNAYGTESFRSYETADRATFEDSAAKFSISRSRTDCTEVSDETTSGI",
               "SFKTEPFGPPSSPESTSDSKITRNLDDCSGCGRQIQDRFYLSAVEKRWHASCLQCYACRQ",
               "PLERESSCYSRDGNIYCKNDYYSFFGTRRCSRCLASISSNELVMRARNLVFHVNCFCCTV",
               "CHTPLTKGDQYGIIDALIYCRTHYSIAREGDTASSSMSATYPYSAQFGSPHNDSSSPHSD",
               "PSRSIVPTGIFVPASHVINGLPQPARQKGRPRKRKPKDIEAFTANIDLNTEYVDFGRGSH",
               "LSSSSRTKRMRTSFKHHQLRTMKSYFAINHNPDAKDLKQLSQKTGLPKRVLQVWFQNARA",
               "KWRRMMMKQDGSGLLEKGEGALDLDSISVHSPTSFILGGPNSTPPLNLD"]
        reviewed = "Reviewed"
        family = "Homeodomain"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # FBgn0000054
        gene = "Adf1"
        uniacc = "P05552"
        unientry = "ADF1_DROME"
        seq = ["MHTLTAAIEMDKLDANLEQQFDLNLIEAVKLNPVIYDRSHYNYKHFVRKAQTWKQIAETL",
               "GVPEQKCTKRWKSLRDKFAREMKLCQESRWRYFKQMQFLVDSIRQYRESLLGKCANGSQS",
               "ANQVADPSQQQQAQQQTVVDIFAQPFNGSATTSAQALTHPHEITVTSDAQLATAVGKDQK",
               "PYFYEPPLKRERSEEEHSDNMLNTIKIFQNNVSQAVSAEDQSFGMVVTDMLNTLGVRQKA",
               "EAKVHIIKYLTDMQLLAQHNKY"]
        reviewed = "Reviewed"
        family = "MADF"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # FBgn0037645
        gene = "ich"
        uniacc = "Q9VHJ6"
        unientry = "ICH_DROME"
        seq = ["MKFCSMNGPSESEVESLLMSPCWGPPQTSAQDQYLLDGHKLLLEHDLDSLPSDADQKNSL",
               "DVLDNLLLNGSSLNALSDLKPLPPFTGYTGHLSINGISGHHYHAIAQRLPDESNNNYMQS",
               "AYHPQNQSNPTSTTQSNGGSNSNSNNSNEHNIVSSSTCLPESVLSGSGGGGGGNDACLAD",
               "VKLFADSQLDSKLYSVADSCVMNGSGSGSAANGNGSGPSIATLTPIDQVADSLHNSGVRI",
               "YKDLTEYVDMNSIDDIAAIIGSAIADTTVPNQLDKDDNNDTRDSWMDLDAWIDGNCIQQE",
               "SAKVLVSQQDSLGDFILPHSPLPMHASSSTLQSLLSHGYMPLLQNRLQNGPPNGNSSGGG",
               "GGANQGAGIKGDPQAPTSTSYCNELAAATSSSCSPPGSVVSTTDNPNGLMINPRYLSNPT",
               "NNQATGNLHQQGGYSMQASGLSLKDNGLCSPDLLGNYPHTTTASTTGSEMRTGAPKAKRS",
               "RSQKKSNQQQQQQQQQQQQQGDGGGQPTTPQMSAISPSGFSASDLSGLLGKEKPVHRCSI",
               "CNRGFLNKSNIKVHLRTHTGEKPFRCDVCAKAFRQKAHLLKHQQIHKRIGRD"]
        reviewed = "Reviewed"
        family = "C2H2 ZF"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
    elif s == "Homo_sapiens":
        # ENSG00000214189
        gene = "ZNF788P"
        uniacc = "Q6ZQV5"
        unientry = "ZN788_HUMAN"
        seq = ["MRNMIPQDNENPPQQGEANQNDSVAFEDVAVNFTPDEWALLDPSQKNLYREVMQETLRNL",
               "ASIEVLWKRDSLKVKVISMEKF"]
        reviewed = "Reviewed"
        family = "C2H2 ZF"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # ENSG00000249459
        gene = "ZNF286B"
        uniacc = "P0CG31"
        unientry = "Z286B_HUMAN"
        seq = ["METDLAEMPEKGVLSSQDSPHFQEKSTEEGEVAALRLTARSQAAAAAAAPGSRSLRGVHV",
               "PPPLHPAPAREEIKSTCSLKACFSLSLTLTYYRTAFLLSTENEGNLHFQCPSDVETRPQS",
               "KDSTSVQDFSKAESCKVAIIDRLTRNSVYDSNLEAALECENWLEKQQGNQERHLREMFTH",
               "MNSLSEETDHEHDVYWKSFNQKSVLITEDRVPKGSYAFHTLEKSLKQKSNLMKKQRTYKE",
               "KKPHKCNDCGELFTCHSVHIQHQRVHTGEKPYTCNECGKSFSHRANLTKHQRTHTRILFE",
               "CRECKKTFTESSSLATHQRIHVGERPYECNECGKGFNRSTHLVQHQLIHTGVRPYECNEC",
               "DKAFIHSSALIKHQRTHTGEKPYKCQECGKAFSHCSSLTKHQRVHTGEKPYECSECGKTF",
               "SQSTHLVQHQRIHTGEKPYECSECGKTFSQSSNFAKHQRIHIGKKPYKCSECGKAFIHSS",
               "ALIQHQRTHTGEKPFRCNECGKSFKCSSSLIRHQRVHTEEQP"]
        reviewed = "Reviewed"
        family = "C2H2 ZF"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # DUX1_HUMAN
        gene = "DUX1"
        uniacc = "O43812"
        unientry = "DUX1_HUMAN"
        seq = ["MALLTALDDTLPEEAQGPGRRMILLSTPSQSDALRACFERNLYPGIATKEELAQGIDIPE",
               "PRVQIWFQNERSCQLRQHRRQSRPWPGRRDPQKGRRKRTAITGSQTALLLRAFEKDRFPG",
               "IAAREELARETGLPESRIQIWFQNRRARHRGQSGRAPTQASIRCNAAPIG"]
        reviewed = "Reviewed"
        family = "Homeodomain"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # DUX3_HUMAN
        gene = "DUX3"
        uniacc = "Q96PT4"
        unientry = "DUX3_HUMAN"
        seq = ["MPAEVHGSPPASLCPCPSVKFRPGLPAMALLTALDDTLPEEAQGPGRRMILLSTPSQSDA",
               "LRACFERNLYPGIATKEQLAQGIDIPEPRVQIWFQNERSCQLRQHRRQSRPWPGRRDPQK",
               "GRRKRTAITGSQTALLLRAFEKDRFPGIPAREELARETGLPESRIQLWFQNRRARHWGQS",
               "GRAPTQASIRCNAAPIG"]
        reviewed = "Reviewed"
        family = "Homeodomain"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # ENSG00000250709
        gene = "SOHLH2"
        uniacc = "Q9NX45"
        unientry = "SOLH2_HUMAN"
        seq = ["MASSIICQEHCQISGQAKIDILLVGDVTVGYLADTVQKLFANIAEVTITISDTKEAAALL",
               "DDCIFNMVLLKVPSSLSAEELEAIKLIRFGKKKNTHSLFVFIIPENFKGCISGHGMDIAL",
               "TEPLTMEKMSNVVKYWTTCPSNTVKTENATGPEELGLPLQRSYSEHLGYFPTDLFACSES",
               "LRNGNGLELNASLSEFEKNKKISLLHSSKEKLRRERIKYCCEQLRTLLPYVKGRKNDAAS",
               "VLEATVDYVKYIREKISPAVMAQITEALQSNMRFCKKQQTPIELSLPGTVMAQRENSVMS",
               "TYSPERGLQFLTNTCWNGCSTPDAESSLDEAVRVPSSSASENAIGDPYKTHISSAALSLN",
               "SLHTVRYYSKVTPSYDATAVTNQNISIHLPSAMPPVSKLLPRHCTSGLGQTCTTHPNCLQ",
               "QFWAY"]
        reviewed = "Reviewed"
        family = "bHLH"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # ENSG00000233608
        gene = "TWIST2"
        uniacc = "Q8WVJ9"
        unientry = "TWST2_HUMAN"
        seq = ["MEEGSSSPVSPVDSLGTSEEELERQPKRFGRKRRYSKKSSEDGSPTPGKRGKKGSPSAQS",
               "FEELQSQRILANVRERQRTQSLNEAFAALRKIIPTLPSDKLSKIQTLKLAARYIDFLYQV",
               "LQSDEMDNKMTSCSYVAHERLSYAFSVWRMEGAWSMSASH"]
        reviewed = "Reviewed"
        family = "bHLH"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
    elif s == "Mus_musculus":
        # BAE29047.1
        gene = "Cebpd"
        uniacc = "Q3UE64"
        unientry = "Q3UE64_MOUSE"
        seq = ["SPGPSLAPGTVREKGAGKRGPDRGSPEYRQRRERNNIAVRKSRDKAKRRNQEMQQKLVEL",
               "SAENEKLHQRVEQLTRDLAGLRQFFKKLPSPPFLPPTGADCR"]
        reviewed = "Unreviewed"
        family = "bZIP"
        jaspar_motifs = ""
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # A1JVI6_MOUSE
        gene = "Dux4"
        uniacc = "A1JVI6"
        unientry = "A1JVI6_MOUSE"
        seq = ["MAEAGSPVGGSGVARESRRRRKTVWQAWQEQALLSTFKKKRYLSFKERKELAKRMGVSDC",
               "RIRVWFQNRRNRSGEEGHASKRSIRGSRRLASPQLQEELGSRPQGRGMRSSGRRPRTRLT",
               "SLQLRILGQAFERNPRPGFATREELARDTGLPEDTIHIWFQNRRARRRHRRGRPTAQDQD",
               "LLASQGSDGAPAGPEGREREGAQENLLPQEEAGSTGMDTSSPSDLPSFCGESQPFQVAQP",
               "RGAGQKEDPTRAGNAGSLEPLLDQLLDEVQVEEPAPAPLNLDGDPGGRVHEGSQESFWPQ",
               "EEAGSTGMDTSSPSDSNSFCRESQPSQVAQPCGAGQEDARTQADSTGPLELLLLDQLLDE",
               "VQKEEHVPVPLDWGRNPGSREHEGSQDSLLPLEEAVNSGMDTSIPSIWPTFCRESQPPQV",
               "AQPSGPGQAQAPTQGGNTDPLELFLYQLLDEVQVEEHAPAPLNWDVDPGGRVHEGSWESF",
               "WPQEEAGSTGLDTSSPSDSNSFFRESKPSQVAQRRGAGQEDARTQADSTGPLELLLFDQL",
               "LDEVQKEEHVPAPLDWGRNPGSMEHEGSQDSLLPLEEAANSGRDTSIPSIWPAFCRKSQP",
               "PQVAQPSGPGQAQAPIQGGNTDPLELFLDQLLTEVQLEEQGPAPVNVEETWEQMDTTPDL",
               "PLTSEEYQTLLDML"]
        reviewed = "Unreviewed"
        family = "Homeodomain"
        jaspar_motifs = ""
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # ENSMUSG00000073448
        gene = "mszf57"
        uniacc = "O88255"
        unientry = "O88255_MOUSE"
        seq = ["EKPYKCKECGKSFRQTSVLKSHQKMHTGEKPYKCKQCDKSFAHSSSFRTHQKIHTSEEHC",
               "SCPECGREFHQLSHLRKHYRLHTGE"]
        reviewed = "Unreviewed"
        family = "C2H2 ZF"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))     
        # Q8K439_MOUSE
        gene = "Zfp263"
        uniacc = "Q8K439"
        unientry = "Q8K439_MOUSE"
        seq = ["MTMAAGPSSQEPEGLLIVKLEEDCAWSHEVPPPEPEPSPEASHLRFRRFRFQDAPGPREA",
               "LSRLQELCRGWLRPEMRTKEQILELLVLEQFLTILPQEIQSRVQELRPESGEEAVILVER",
               "MQKELGKLRQQFTNQGRGAEVLLEEPLPLETAGESPSFKLEPMEIERSPGPRLQELLDPS",
               "PQRDSQAVKERALSAPWLSLFPPEGNVEDKDMTGTQLPESLEDMAMYISQEWDHQDPSKR",
               "ALSRYMVQDSYENSGTLESSIPSQEVSSTHVEQGEKLWDSSVQTCKEGMNPRNPVPGVEK",
               "FENQERNVESVSPESTHPPVLLPGQARREVPWSPEQGRLDDREGHWECPPEDKIEESLVG",
               "TPSCKGLVQAKEQPKKLHLCALCGKNFSNNSNLIRHQRIHAAEKLCMDVECGEVFGGHPH",
               "FLSLHRTHIGEEAHKCLECGKCFSQNTHLTRHQRTHTGEKPFQCNACGKSFSCNSNLNRH",
               "QRTHTGEKPYKCPECGEIFAHSSNLLRHQRIHTGERPYRCSECGKSFSRSSHLVIHERTH",
               "EKERLDPFPECGQGMNDSAPFLTNHRVEKKLFECSTCGKSFRQGMHLTRHQRTHTGEKPY",
               "KCILCGENFSHRSNLIRHQRIHTGEKPYTCHECGDSFSHSSNRIRHLRTHTGERPYKCSE",
               "CGESFSRSSRLTSHQRTHTG"]
        reviewed = "Unreviewed"
        family = "C2H2 ZF"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
    elif s == "Saccharomyces_cerevisiae":
        # MAL63_YEASX
        gene = "MAL63"
        uniacc = "P10508"
        unientry = "MAL63_YEASX"
        seq = ["MGIAKQSCDCCRVRRVKCDRNKPCNRCIQRNLNCTYLQPLKKRGPKSIRAGSLKKIAEVQ",
               "MVSMNNNIMAAPVVCKKVPKNLIDQCLRLRLYHDNLYVIWPMLSYDDLHKLLEEKYDDRC",
               "AYWFLVSLSAATLSDLQIEIEYEEGVTFTGEQLCTLCMLSRQFFDDLSNSDIFRIMTYYC",
               "LHRCYAQFADTRTSYRLSCEAVGLIIKIAGFHREETYEFLPFGEQQLRRKVYYLLLMTER",
               "FYAVYIKCVTSLDATIAPPLPEVVTDPRLSLESFLEVIRVFTIPGKCFYDALATNCVDDS",
               "CTEDSLKRIRNELHTTSLDIEPWSYGYIDFLFSRHWVRTLAWKLVLHMKGMRMNFLSNTN",
               "NTHIPVEIARDMLGDTFLTPKNLYDVHGPGIPMKALEIANALVDVVNKYDHNMKLEAWNV",
               "LYDVSKFVFSLKHCNNKMFDRFSTKCQGALITLPISKPLQLNDNSKDEDDIIP"]
        reviewed = "Reviewed"
        family = "Zinc cluster"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # YML007W
        gene = "YAP1"
        uniacc = "P19880"
        unientry = "AP1_YEAST"
        seq = ["MSVSTAKRSLDVVSPGSLAEFEGSKSRHDEIENEHRRTGTRDGEDSEQPKKKGSKTSKKQ",
               "DLDPETKQKRTAQNRAAQRAFRERKERKMKELEKKVQSLESIQQQNEVEATFLRDQLITL",
               "VNELKKYRPETRNDSKVLEYLARRDPNLHFSKNNVNHSNSEPIDTPNDDIQENVKQKMNF",
               "TFQYPLDNDNDNDNSKNVGKQLPSPNDPSHSAPMPINQTQKKLSDATDSSSATLDSLSNS",
               "NDVLNNTPNSSTSMDWLDNVIYTNRFVSGDDGSNSKTKNLDSNMFSNDFNFENQFDEQVS",
               "EFCSKMNQVCGTRQCPIPKKPISALDKEVFASSSILSSNSPALTNTWESHSNITDNTPAN",
               "VIATDATKYENSFSGFGRLGFDMSANHYVVNDNSTGSTDSTGSTGNKNKKNNNNSDDVLP",
               "FISESPFDMNQVTNFFSPGSTGIGNNAASNTNPSLLQSSKEDIPFINANLAFPDDNSTNI",
               "QLQPFSESQSQNKFDYDMFFRDSSKEGNNLFGEFLEDDDDDKKAANMSDDESSLIKNQLI",
               "NEEPELPKQYLQSVPGNESEISQKNGSSLQNADKINNGNDNDNDNDVVPSKEGSLLRCSE",
               "IWDRITTHPKYSDIDVDGLCSELMAKAKCSERGVVINAEDVQLALNKHMN"]
        reviewed = "Reviewed"
        family = "bZIP"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))
        # YDR423C
        gene = "CAD1"
        uniacc = "P24813"
        unientry = "AP2_YEAST"
        seq = ["MGNILRKGQQIYLAGDMKKQMLLNKDGTPKRKVGRPGRKRIDSEAKSRRTAQNRAAQRAF",
               "RDRKEAKMKSLQERVELLEQKDAQNKTTTDFLLCSLKSLLSEITKYRAKNSDDERILAFL",
               "DDLQEQQKRENEKGTSTAVSKAAKELPSPNSDENMTVNTSIEVQPHTQENEKVMWNIGSW",
               "NAPSLTNSWDSPPGNRTGAVTIGDESINGSEMPDFSLDLVSNDRQTGLEALDYDIHNYFP",
               "QHSERLTAEKIDTSACQCEIDQKYLPYETEDDTLFPSVLPLAVGSQCNNICNRKCIGTKP",
               "CSNKEIKCDLITSHLLNQKSLASVLPVAASHTKTIRTQSEAIEHISSAISNGKASCYHIL",
               "EEISSLPKYSSLDIDDLCSELIIKAKCTDDCKIVVKARDLQSALVRQLL"]
        reviewed = "Reviewed"
        family = "bZIP"
        jaspar_motifs = set()
        if uniacc in jaspar2uniprot:
            for jaspar_motif in jaspar2uniprot[uniacc][0]:
                jaspar_motifs.add(jaspar_motif)
        lines.add("%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s" % (\
                  gene, ss, uniacc, unientry, "".join(seq), reviewed,
                  family, ";".join(sorted(jaspar_motifs))))

    # For each line...
    for line in sorted(lines):
        print(line)