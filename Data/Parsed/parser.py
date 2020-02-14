import argparse
import coreapi
import json
import os
import re
import sys

# Get root path
root_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        os.pardir, os.pardir)

# Append JASPAR-profile-inference to path
jaspar_dir = os.path.join(root_dir, "JASPAR-profile-inference")
sys.path.insert(0, os.path.realpath(jaspar_dir))

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
        self.family = family
        self.evidence = evidence

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

print("Name\tSpecies\tUniAcc\tUniEntry\tSequence\tStatus\tFamily\tJASPAR")

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
    fas_file = os.path.join(root_dir, "Data", "Databases", "Uniprot",
                            "%s_TFs.fa.gz" % s)
    for seq_record in Jglobals.parse_fasta_file(fas_file):
        m = re.search("^\S+\|(\S+)\|(\S+) .+ OX=\d+ GN=(.+) PE=\d+ ",
                        seq_record.description)
        if m:
            fas.setdefault(m.group(2), [m.group(1), m.group(3),
                                        str(seq_record.seq)])

    # For each line...
    txt = {}
    txt_file = os.path.join(root_dir, "Data", "Databases", "Uniprot",
                            "%s_TFs.txt.gz" % s)
    for line in Jglobals.parse_file(txt_file):
        if line.startswith("//"):
            if unientry in fas:
                for gene in genes:
                    txt.setdefault(gene, {})
                    txt[gene].setdefault(reviewed, [])
                    # txt[gene][reviewed].append((uniaccs, unientry, orthodb))
                    txt[gene][reviewed].append((uniaccs, unientry))
        if line.startswith("ID"):
            m = re.search("^ID\s+(\S+)\s+(Reviewed|Unreviewed);", line)
            unientry = m.group(1)
            reviewed = m.group(2)
            uniaccs = []
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
    lines = []

    # For each CIS-BP object...
    for cisbp in sorted(tfs):

        if cisbp.geneid not in txt:
            continue

        # Initialize
        uniaccs = []
        unientries = []
        family = cisbp.family
        sequences = []
        jaspar_motifs = set()
        for reviewed in ["Reviewed", "Unreviewed"]:
            if reviewed not in txt[cisbp.geneid]:
                continue
            for uaccs, uentry in txt[cisbp.geneid][reviewed]:
                for uniacc in uaccs:
                    if uniacc in jaspar2uniprot:
                        for jaspar_motif in jaspar2uniprot[uniacc][0]:
                            jaspar_motifs.add(jaspar_motif)
                if fas[uentry][0] not in uniaccs:
                    uniaccs.append(fas[uentry][0])
                    unientries.append(uentry)
                    sequences.append(fas[uentry][2])
            break
        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (\
                     fas[uentry][1], ss, ";".join(uniaccs), ";".join(unientries),
                     ";".join(sequences), reviewed, family, ";".join(sorted(jaspar_motifs))))

    # Fix species-specific cases
    if s == "Arabidopsis_thaliana":
        # T4P13.29
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
        jaspar_motifs = ""
        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (\
                     gene, ss, uniacc, unientry, ";".join(seq), reviewed,
                     family, jaspar_motifs))
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
        jaspar_motifs = ""
        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (\
                     gene, ss, uniacc, unientry, ";".join(seq), reviewed,
                     family, jaspar_motifs))
    elif s == "Caenorhabditis_elegans":
        # hlh-27
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
        jaspar_motifs = ""
        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (\
                     gene, ss, uniacc, unientry, ";".join(seq), reviewed,
                     family, jaspar_motifs))
    elif s == "Drosophila_melanogaster":
        # br
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
        jaspar_motifs = "MA0010.1;MA0011.1;MA0012.1"
        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (\
                     gene, ss, uniacc, unientry, ";".join(seq), reviewed,
                     family, jaspar_motifs))
    elif s == "Homo_sapiens":
        # DUX1
        gene = "DUX1"
        uniacc = "O43812"
        unientry = "DUX1_HUMAN"
        seq = ["MALLTALDDTLPEEAQGPGRRMILLSTPSQSDALRACFERNLYPGIATKEELAQGIDIPE",
               "PRVQIWFQNERSCQLRQHRRQSRPWPGRRDPQKGRRKRTAITGSQTALLLRAFEKDRFPG",
               "IAAREELARETGLPESRIQIWFQNRRARHRGQSGRAPTQASIRCNAAPIG"]
        reviewed = "Reviewed"
        family = "Homeodomain"
        jaspar_motifs = ""
        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (\
                     gene, ss, uniacc, unientry, ";".join(seq), reviewed,
                     family, jaspar_motifs))
        # DUX3
        gene = "DUX3"
        uniacc = "Q96PT4"
        unientry = "DUX3_HUMAN"
        seq = ["MPAEVHGSPPASLCPCPSVKFRPGLPAMALLTALDDTLPEEAQGPGRRMILLSTPSQSDA",
               "LRACFERNLYPGIATKEQLAQGIDIPEPRVQIWFQNERSCQLRQHRRQSRPWPGRRDPQK",
               "GRRKRTAITGSQTALLLRAFEKDRFPGIPAREELARETGLPESRIQLWFQNRRARHWGQS",
               "GRAPTQASIRCNAAPIG"]
        reviewed = "Reviewed"
        family = "Homeodomain"
        jaspar_motifs = ""
        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (\
                     gene, ss, uniacc, unientry, ";".join(seq), reviewed,
                     family, jaspar_motifs))
    elif s == "Mus_musculus":
        # Dux4
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
        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (\
                     gene, ss, uniacc, unientry, ";".join(seq), reviewed,
                     family, jaspar_motifs))
        # Cebpd
        gene = "Cebpd"
        uniacc = "Q3UE64"
        unientry = "Q3UE64_MOUSE"
        seq = ["SPGPSLAPGTVREKGAGKRGPDRGSPEYRQRRERNNIAVRKSRDKAKRRNQEMQQKLVEL",
               "SAENEKLHQRVEQLTRDLAGLRQFFKKLPSPPFLPPTGADCR"]
        reviewed = "Unreviewed"
        family = "bZIP"
        jaspar_motifs = ""
        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (\
                     gene, ss, uniacc, unientry, ";".join(seq), reviewed,
                     family, jaspar_motifs))
        # Zfp263
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
        jaspar_motifs = ""
        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (\
                     gene, ss, uniacc, unientry, ";".join(seq), reviewed,
                     family, jaspar_motifs))
    elif s == "Saccharomyces_cerevisiae":
        # ARG5,6
        gene = "ARG5,6"
        uniacc = "Q01217"
        unientry = "ARG56_YEAST"
        seq = ["MPSASLLVSTKRLNASKFQKFVSSLNKSTIAGFASVPLRAPPSVAFTRKKVGYSKRYVSS",
               "TNGFSATRSTVIQLLNNISTKREVEQYLKYFTSVSQQQFAVIKVGGAIISDNLHELASCL",
               "AFLYHVGLYPIVLHGTGPQVNGRLEAQGIEPDYIDGIRITDEHTMAVVRKCFLEQNLKLV",
               "TALEQLGVRARPITSGVFTADYLDKDKYKLVGNIKSVTKEPIEASIKAGALPILTSLAET",
               "ASGQMLNVNADVAAGELARVFEPLKIVYLNEKGGIINGSTGEKISMINLDEEYDDLMKQS",
               "WVKYGTKLKIREIKELLDYLPRSSSVAIINVQDLQKELFTDSGAGTMIRRGYKLVKRSSI",
               "GEFPSADALRKALQRDAGISSGKESVASYLRYLENSDFVSYADEPLEAVAIVKKDTNVPT",
               "LDKFVCSDAAWLNNVTDNVFNVLRRDFPALQWVVSENDANIAWHFDKSQGSYLKGGKVLF",
               "WYGIDDINTISELVENFVKSCDTASTLNSSASSGVFANKKSARSYSTRSTPRPEGVNTNP",
               "GRVALIGARGYTGKNLVSLINGHPYLEVAHVSSRELKGQKLQDYTKSEIIYESLQIQDIR",
               "KLEEQNAVDFWVMALPNKVCEPFVETIQSVHGKSKIIDLSADHRFVSESDWAYGLPELND",
               "RAKIANAAKIANPGCYATGSQLTISPLTKYINGLPTVFGVSGYSGAGTKPSPKNDPKFLN",
               "NNLIPYALSDHIHEREISARIGHNVAFMPHVGQWFQGISLTVSIPIKKGSLSIDEIRKLY",
               "RNFYEDEKLVHVIDDIPLVKDIEGTHGVVIGGFKLNDAEDRVVVCATIDNLLKGAATQCL",
               "QNINLAMGYGEYAGIPENKIIGV"]
        reviewed = "Reviewed"
        family = "Unknown"
        jaspar_motifs = ""
        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (\
                     gene, ss, uniacc, unientry, ";".join(seq), reviewed,
                     family, jaspar_motifs))
        # MAL63
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
        jaspar_motifs = ""
        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (\
                     gene, ss, uniacc, unientry, ";".join(seq), reviewed,
                     family, jaspar_motifs))

    # For each line...
    for line in sorted(lines):
        print(line)