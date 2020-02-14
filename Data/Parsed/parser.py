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

#-------------#
# Functions   #
#-------------#

def uniacc_to_OrthoDB_cluster(uniacc):

    # Initialize
    codec = coreapi.codecs.CoreJSONCodec()
    client = coreapi.Client()
    orthodb_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                               ".OrthoDB")

    if not os.path.isdir(orthodb_dir):
        os.mkdir(orthodb_dir)

    # Skip if already JSON data already fetched
    json_file = os.path.join(orthodb_dir, ".%s.json" % uniacc)
    if not os.path.exists(json_file):

        # Get response
        url = "https://www.orthodb.org/search?query=%s&level=2759&species=2759" % \
            uniacc
        response = client.get(url)
        json_obj = json.loads(codec.encode(response))
        with open(json_file, "w") as j:
            j.write(json.dumps(json_obj))

    with open(json_file, "r") as j:  
        json_obj = json.load(j)
        for orthodb_cluster in json_obj["data"]:
            yield(orthodb_cluster)

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

print("Name\tUniAcc\tUniEntry\tSequence\tStatus\tFamily\tOrthoDB\tJASPAR")

for s in sorted(species):

    # Initialize
    lines = []

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
                    txt[gene][reviewed].append((uniaccs, unientry))
        if line.startswith("ID"):
            m = re.search("^ID\s+(\S+)\s+(Reviewed|Unreviewed);", line)
            unientry = m.group(1)
            reviewed = m.group(2)
            uniaccs = []
            genes = set()
            pfams = set()
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
        clusters = set()
        for reviewed in ["Reviewed", "Unreviewed"]:
            if reviewed not in txt[cisbp.geneid]:
                continue
            for uaccs, uentry in txt[cisbp.geneid][reviewed]:
                for uniacc in uaccs:
                    if uniacc in jaspar2uniprot:
                        for jaspar_motif in jaspar2uniprot[uniacc][0]:
                            jaspar_motifs.add(jaspar_motif)
                    for cluster in uniacc_to_OrthoDB_cluster(uniacc):
                        clusters.add(cluster)
                if fas[uentry][0] not in uniaccs:
                    uniaccs.append(fas[uentry][0])
                    unientries.append(uentry)
                    sequences.append(fas[uentry][2])
            break

        lines.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % \
                     (fas[uentry][1], ";".join(uniaccs), ";".join(unientries),
                     ";".join(sequences), reviewed, family, ";".join(sorted(clusters)),
                     ";".join(sorted(jaspar_motifs))))

    # Fix species-specific cases
    if s == "Arabidopsis_thaliana":
        clusters = set()
        for cluster in uniacc_to_OrthoDB_cluster("Q9MAB7"):
            clusters.add(cluster)
        seq = ["MEEQQEFEKPIFEFRPKKLRTSRVSRNLMKKTGFRESMNHYEELSCNYGLRENPKKTQKS",
               "LLDHRFCTRRRRNKKILIRCKECGKGFLYEKCLFNHLQVTHSEESTRRSLFCRFSIVQRR",
               "KRSKRVSRYKKILPRFSVSSSSCTMFPVSVDDGGLLEVAESLILLSMSGGKFVNGLEHFG",
               "KALGSTQRKFEDGLLRNEQRLVGEALDSNPEKKLVGISRASVGTSKELSGYLANKKGRED",
               "DELGQQKQAGARILREETDNEQKLVRQETAFEDSVSGFEMNIEHRCGLCHKVFSTYQTLG",
               "GHQTFHRMRNKSKSQTKRCREESIEAEAGINNGSVTLTISAEFAEGCLGQYML"]
        lines.append("T4P13.29\tQ9MAB7\tQ9MAB7_ARATH\t%s\tUnreviewed\tC2H2 ZF\t%s\t" % \
                     ("".join(seq), ";".join(sorted(clusters))))
        clusters = set()
        for cluster in uniacc_to_OrthoDB_cluster("O49746"):
            clusters.add(cluster)
        seq = ["MGRAPCCDKANVKKGPWSPEEDAKLKSYIENSGTGGNWIALPQKIGLKRCGKSCRLRWLN",
               "YLRPNIKHGGFSEEEENIICSLYLTIGSRWSIIAAQLPGRTDNDIKNYWNTRLKKKLINK",
               "QRKELQEACMEQQEMMVMMKRQHQQQQIQTSFMMRQDQTMFTWPLHHHNVQVPALFRIKP",
               "TRFATKKMLSQCSSRTWSRSKIKNWRKQTSSSSRFNDNAFDHLSFSQLLLDPNHNHLGSG",
               "EGFSMNSILSANTNSPLLNTSNDNQWFGNFQAETVNLFSGASTSTSADQSTISWEDISSL",
               "VYSDSKQFF"]
        lines.append("AtMYB84\tO49746\tO49746_ARATH\t%s\tUnreviewed\tMyb/SANT\t%s\t" % \
                    ("".join(seq), ";".join(sorted(clusters))))
    elif s == "Caenorhabditis_elegans":
        clusters = set()
        for cluster in uniacc_to_OrthoDB_cluster("Q18056"):
            clusters.add(cluster)
        seq = ["MPKVIPSSMSDYRSVPYNQTPKSGSERKRRNITNELINECKTIVQKSEEEHISQEVVLFR",
               "IVKLVTGVNLESNFSSNDLSESTRRKFDTESERRKVKTEREKIRRKKQDDCYAELKFFIL",
               "NKQMGSYEQRLKLERITILEIIIDYIKHNSDLLYPETIPQILPLLAGKSTATCENKENEK",
               "PKTRMEVKDLFPRLTFQEVQESPTSTSPLLTFPCIPMIPTTQFNVLSNYNTVPSIFSAPL",
               "RFILPSLQILTPETSDEEENEETVDISN"]
        lines.append("hlh-27\tQ18056\tQ18056_CAEEL\t%s\tUnreviewed\tbHLH\t%s\t" % \
                     ("".join(seq), ";".join(sorted(clusters))))
    elif s == "Drosophila_melanogaster":
        clusters = set()
        for cluster in uniacc_to_OrthoDB_cluster("Q01295"):
            clusters.add(cluster)
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
        lines.append("br\tQ01295\tBRC1_DROME\t%s\tReviewed\tC2H2 ZF\t%s\tMA0010.1;MA0011.1;MA0012.1" % \
                     ("".join(seq), ";".join(sorted(clusters))))
    elif s == "Homo_sapiens":
        clusters = set()
        for cluster in uniacc_to_OrthoDB_cluster("O43812"):
            clusters.add(cluster)
        seq = ["MALLTALDDTLPEEAQGPGRRMILLSTPSQSDALRACFERNLYPGIATKEELAQGIDIPE",
               "PRVQIWFQNERSCQLRQHRRQSRPWPGRRDPQKGRRKRTAITGSQTALLLRAFEKDRFPG",
               "IAAREELARETGLPESRIQIWFQNRRARHRGQSGRAPTQASIRCNAAPIG"]
        lines.append("DUX1\tO43812\tDUX1_HUMAN\t%s\tReviewed\tHomeodomain\t%s\t" % \
                     ("".join(seq), ";".join(sorted(clusters))))
        clusters = set()
        for cluster in uniacc_to_OrthoDB_cluster("Q96PT4"):
            clusters.add(cluster)
        seq = ["MPAEVHGSPPASLCPCPSVKFRPGLPAMALLTALDDTLPEEAQGPGRRMILLSTPSQSDA",
               "LRACFERNLYPGIATKEQLAQGIDIPEPRVQIWFQNERSCQLRQHRRQSRPWPGRRDPQK",
               "GRRKRTAITGSQTALLLRAFEKDRFPGIPAREELARETGLPESRIQLWFQNRRARHWGQS",
               "GRAPTQASIRCNAAPIG"]
        lines.append("DUX3\tQ96PT4\tDUX3_HUMAN\t%s\tReviewed\tHomeodomain\t%s\t" % \
                     ("".join(seq), ";".join(sorted(clusters))))
    elif s == "Mus_musculus":
        clusters = set()
        for cluster in uniacc_to_OrthoDB_cluster("A1JVI6"):
            clusters.add(cluster)
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
        lines.append("Dux4\tA1JVI6\tA1JVI6_MOUSE\t%s\tUnreviewed\tHomeodomain\t%s\t" % \
                     ("".join(seq), ";".join(sorted(clusters))))
        clusters = set()
        for cluster in uniacc_to_OrthoDB_cluster("Q3UE64"):
            clusters.add(cluster)
        seq = ["SPGPSLAPGTVREKGAGKRGPDRGSPEYRQRRERNNIAVRKSRDKAKRRNQEMQQKLVEL",
               "SAENEKLHQRVEQLTRDLAGLRQFFKKLPSPPFLPPTGADCRL"]
        lines.append("Cebpd\tQ3UE64\tQ3UE64_MOUSE\t%s\tUnreviewed\tbZIP\t%s\t" % \
                     ("".join(seq), ";".join(sorted(clusters))))
        clusters = set()
        for cluster in uniacc_to_OrthoDB_cluster("Q8K439"):
            clusters.add(cluster)
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
        lines.append("Zfp263\tQ8K439\tQ8K439_MOUSE\t%s\tUnreviewed\tC2H2 ZF\t%s\t" % \
                     ("".join(seq), ";".join(sorted(clusters))))
    elif s == "Saccharomyces_cerevisiae":
        clusters = set()
        for cluster in uniacc_to_OrthoDB_cluster("Q01217"):
            clusters.add(cluster)
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
        lines.append("ARG5,6\tQ01217\tARG56_YEAST\t%s\tReviewed\tUnknown\t%s\t" % \
                     ("".join(seq), ";".join(sorted(clusters))))
        clusters = set()
        for cluster in uniacc_to_OrthoDB_cluster("P10508"):
            clusters.add(cluster)
        seq = ["MGIAKQSCDCCRVRRVKCDRNKPCNRCIQRNLNCTYLQPLKKRGPKSIRAGSLKKIAEVQ",
               "MVSMNNNIMAAPVVCKKVPKNLIDQCLRLRLYHDNLYVIWPMLSYDDLHKLLEEKYDDRC",
               "AYWFLVSLSAATLSDLQIEIEYEEGVTFTGEQLCTLCMLSRQFFDDLSNSDIFRIMTYYC",
               "LHRCYAQFADTRTSYRLSCEAVGLIIKIAGFHREETYEFLPFGEQQLRRKVYYLLLMTER",
               "FYAVYIKCVTSLDATIAPPLPEVVTDPRLSLESFLEVIRVFTIPGKCFYDALATNCVDDS",
               "CTEDSLKRIRNELHTTSLDIEPWSYGYIDFLFSRHWVRTLAWKLVLHMKGMRMNFLSNTN",
               "NTHIPVEIARDMLGDTFLTPKNLYDVHGPGIPMKALEIANALVDVVNKYDHNMKLEAWNV",
               "LYDVSKFVFSLKHCNNKMFDRFSTKCQGALITLPISKPLQLNDNSKDEDDIIP"]
        lines.append("MAL63\tP10508\tMAL63_YEASX\t%s\tReviewed\tZinc cluster\t%s\t" % \
                     ("".join(seq), ";".join(sorted(clusters))))

    # For each line...
    for line in sorted(lines):
        print(line)