import coreapi
import json

client = coreapi.Client()
codec = coreapi.codecs.CoreJSONCodec()

print("tf\tspecies\tdataset")

for taxid in [9606, 3702]:

    tfs = {}
    url = "http://remap.univ-amu.fr:80/api/v1/datasets/findByTaxid/taxid=%s" % taxid
    response = client.get(url)
    remap = json.loads(codec.encode(response))
    for dataset in remap["datasets"]:
        tf = dataset["target_name"]
        tfs.setdefault(tf, [])
        tfs[tf].append(dataset["dataset_name"])
    for tf in sorted(tfs):
        print("%s\t%s\t%s" % (tf, remap["species"]["tax_name"], ";".join(tfs[tf])))