#!/usr/bin/env python

import sys
import csv
import hashlib
import os
from os.path import expandvars
import pprint
from Levenshtein import distance, editops

known_references = {
    "027f4cb7c33de171c14c83c69b90eea6": "hrv-a1 (NC_038311.1)",
    "105c82802b67521950854a851fc6eefd": "sars-cov-2 (MN908947.3)",
    "1b0da1dbf9374d31d98f26e9edbdfe86": "hrv-b3 (NC_038312.1)",
    "2cd494e3006b363034acb3988a135b4b": "hcov-229e (AF304460.1)",
    "33db59ad4cba6c4be72da7bf64b81444": "hrv-c (NC_009996.1)",
    "44be42f112c97e6a84d690636e28b13f": "hcov-nl63 (AY567487.2)",
    "4a7cb2305253cf8809bde534ae007c71": "hrv-b14 (NC_001490.1)",
    "9cf54db346dedfb6b43dfd68b4c7bc51": "hcov-hku1 (AY597011.2)",
    "b9a2b7c3ce1194ebc9ffadca232f548f": "hrv-a (NC_001617.1)",
    "f34f1e552c7a5e94ac72c509b93ac2c4": "hcov-oc43 (AY585228.1)",
}

reference_to_path = {
    "027f4cb7c33de171c14c83c69b90eea6": "${PIPELINE_ROOT}/reference/hrv-a1/hrv-a1.fna",
    "105c82802b67521950854a851fc6eefd": "${PIPELINE_ROOT}/reference/sars-cov-2/GCA_009858895.3_ASM985889v3_genomic.fna",
    "1b0da1dbf9374d31d98f26e9edbdfe86": "${PIPELINE_ROOT}/reference/hrv-b3/hrv-b3.fna",
    "2cd494e3006b363034acb3988a135b4b": "${PIPELINE_ROOT}/reference/hcov-229e/GCA_000853505.1_ViralProj14913_genomic.fna",
    "33db59ad4cba6c4be72da7bf64b81444": "${PIPELINE_ROOT}/reference/hrv-c/hrv-c.fna",
    "44be42f112c97e6a84d690636e28b13f": "${PIPELINE_ROOT}/reference/hcov-nl63/GCA_000853865.1_ViralProj14960_genomic.fna",
    "4a7cb2305253cf8809bde534ae007c71": "${PIPELINE_ROOT}/reference/hrv-b14/hrv-b14.fna",
    "9cf54db346dedfb6b43dfd68b4c7bc51": "${PIPELINE_ROOT}/reference/hcov-hku1/GCA_000858765.1_ViralProj15139_genomic.fna",
    "b9a2b7c3ce1194ebc9ffadca232f548f": "${PIPELINE_ROOT}/reference/hrv-a/hrv-a.fna",
    "f34f1e552c7a5e94ac72c509b93ac2c4": "${PIPELINE_ROOT}/reference/hcov-oc43/GCA_003972325.1_ASM397232v1_genomic.fna",
}

supported_references = {
    "AF304460.1": "2cd494e3006b363034acb3988a135b4b",
    "AY567487.2": "44be42f112c97e6a84d690636e28b13f",
    "AY585228.1": "f34f1e552c7a5e94ac72c509b93ac2c4",
    "AY597011.2": "9cf54db346dedfb6b43dfd68b4c7bc51",
    "MN908947.3": "105c82802b67521950854a851fc6eefd",
    "NC_001490.1": "4a7cb2305253cf8809bde534ae007c71",
    "NC_001617.1": "b9a2b7c3ce1194ebc9ffadca232f548f",
    "NC_009996.1": "33db59ad4cba6c4be72da7bf64b81444",
    "NC_038311.1": "027f4cb7c33de171c14c83c69b90eea6",
    "NC_038312.1": "1b0da1dbf9374d31d98f26e9edbdfe86",
    "hcov-229e": "2cd494e3006b363034acb3988a135b4b",
    "hcov-229e (AF304460.1)": "2cd494e3006b363034acb3988a135b4b",
    "hcov-hku1": "9cf54db346dedfb6b43dfd68b4c7bc51",
    "hcov-hku1 (AY597011.2)": "9cf54db346dedfb6b43dfd68b4c7bc51",
    "hcov-nl63": "44be42f112c97e6a84d690636e28b13f",
    "hcov-nl63 (AY567487.2)": "44be42f112c97e6a84d690636e28b13f",
    "hcov-oc43": "f34f1e552c7a5e94ac72c509b93ac2c4",
    "hcov-oc43 (AY585228.1)": "f34f1e552c7a5e94ac72c509b93ac2c4",
    "hrv-a": "b9a2b7c3ce1194ebc9ffadca232f548f",
    "hrv-a (NC_001617.1)": "b9a2b7c3ce1194ebc9ffadca232f548f",
    "hrv-a1": "027f4cb7c33de171c14c83c69b90eea6",
    "hrv-a1 (NC_038311.1)": "027f4cb7c33de171c14c83c69b90eea6",
    "hrv-b14": "4a7cb2305253cf8809bde534ae007c71",
    "hrv-b14 (NC_001490.1)": "4a7cb2305253cf8809bde534ae007c71",
    "hrv-b3": "1b0da1dbf9374d31d98f26e9edbdfe86",
    "hrv-b3 (NC_038312.1)": "1b0da1dbf9374d31d98f26e9edbdfe86",
    "hrv-c": "33db59ad4cba6c4be72da7bf64b81444",
    "hrv-c (NC_009996.1)": "33db59ad4cba6c4be72da7bf64b81444",
    "sars-cov-2": "105c82802b67521950854a851fc6eefd",
    "sars-cov-2 (MN908947.3)": "105c82802b67521950854a851fc6eefd",
}

reference_to_organism_accession = {
    "AF304460.1": ["hcov-229e", "AF304460.1"],
    "AY567487.2": ["hcov-nl63", "AY567487.2"],
    "AY585228.1": ["hcov-oc43", "AY585228.1"],
    "AY597011.2": ["hcov-hku1", "AY597011.2"],
    "MN908947.3": ["sars-cov-2", "MN908947.3"],
    "NC_001490.1": ["hrv-b14", "NC_001490.1"],
    "NC_001617.1": ["hrv-a", "NC_001617.1"],
    "NC_009996.1": ["hrv-c", "NC_009996.1"],
    "NC_038311.1": ["hrv-a1", "NC_038311.1"],
    "NC_038312.1": ["hrv-b3", "NC_038312.1"],
    "hcov-229e": ["hcov-229e", "AF304460.1"],
    "hcov-229e (AF304460.1)": ["hcov-229e", "AF304460.1"],
    "hcov-hku1": ["hcov-hku1", "AY597011.2"],
    "hcov-hku1 (AY597011.2)": ["hcov-hku1", "AY597011.2"],
    "hcov-nl63": ["hcov-nl63", "AY567487.2"],
    "hcov-nl63 (AY567487.2)": ["hcov-nl63", "AY567487.2"],
    "hcov-oc43": ["hcov-oc43", "AY585228.1"],
    "hcov-oc43 (AY585228.1)": ["hcov-oc43", "AY585228.1"],
    "hrv-a": ["hrv-a", "NC_001617.1"],
    "hrv-a (NC_001617.1)": ["hrv-a", "NC_001617.1"],
    "hrv-a1": ["hrv-a1", "NC_038311.1"],
    "hrv-a1 (NC_038311.1)": ["hrv-a1", "NC_038311.1"],
    "hrv-b14": ["hrv-b14", "NC_001490.1"],
    "hrv-b14 (NC_001490.1)": ["hrv-b14", "NC_001490.1"],
    "hrv-b3": ["hrv-b3", "NC_038312.1"],
    "hrv-b3 (NC_038312.1)": ["hrv-b3", "NC_038312.1"],
    "hrv-c": ["hrv-c", "NC_009996.1"],
    "hrv-c (NC_009996.1)": ["hrv-c", "NC_009996.1"],
    "sars-cov-2": ["sars-cov-2", "MN908947.3"],
    "sars-cov-2 (MN908947.3)": ["sars-cov-2", "MN908947.3"],
}

results = []


def read_fasta(fa):
    with open(fa, "r") as f:
        lines = []
        name = ""

        while line := f.readline().rstrip():
            if line.startswith(">"):
                name = line.removeprefix(">")
            else:
                lines.append(line)

        fasta = "".join(lines)
        return {
            "name": name,
            "hash": hashlib.md5(fasta.encode("utf-8")).hexdigest(),
            "fasta": fasta,
            "length": len(fasta),
        }


def usage():
    print("usage: " + sys.argv[0] + " <reference> <fa> [<fa>...]")
    print("example:\n\t" + sys.argv[0] + "'hcov-hku1 (AY597011.2)' s1-hku1.fa s2-hku1.fa")
    print("\nor\n")
    print("usage: " + sys.argv[0] + " <reference>")
    print("example:\n\t" + sys.argv[0] + "'hcov-hku1 (AY597011.2)'")

    print("references supported:\n\t" + ", ".join(sorted(supported_references.keys())))
    exit(1)


#
# ./make-reference.py ${PIPELINE_ROOT} $(find ${PIPELINE_ROOT}/reference -name '*.fna' | grep -v panel | grep -iP 'sars|hcov|hrv' | sort)
#
if sys.argv[0].endswith("make-reference.py"):
    hash_to_ref = {}
    hash_to_path = {}
    ref_to_hash = {}
    ref_to_org_acc = {}

    replacement = sys.argv[1]

    for fa in range(2, len(sys.argv)):
        fasta = read_fasta(sys.argv[fa])

        org = os.path.dirname(sys.argv[fa]).split("/")[-1]
        ref = fasta["name"].split(" ")[0]
        key = org + " (" + ref + ")"

        hash_to_path[fasta["hash"]] = os.path.realpath(sys.argv[fa]).replace(replacement, "${PIPELINE_ROOT}")
        hash_to_ref[fasta["hash"]] = key

        ref_to_hash[key] = fasta["hash"]
        ref_to_hash[org] = fasta["hash"]
        ref_to_hash[ref] = fasta["hash"]

        ref_to_org_acc[key] = [org, ref]
        ref_to_org_acc[org] = [org, ref]
        ref_to_org_acc[ref] = [org, ref]

    print("known_references = ", end="")
    pprint.pprint(hash_to_ref)
    print("")

    print("reference_to_path = ", end="")
    pprint.pprint(hash_to_path)
    print("")

    print("supported_references = ", end="")
    pprint.pprint(ref_to_hash)
    print("")

    print("reference_to_organism_accession = ", end="")
    pprint.pprint(ref_to_org_acc)
    print("")

    exit(0)


if len(sys.argv) == 2:
    if sys.argv[1] in supported_references:
        print(supported_references[sys.argv[1]])
    else:
        usage()

elif len(sys.argv) < 3:
    usage()

target = supported_references[sys.argv[1]]
target_fasta = read_fasta(expandvars(reference_to_path[target]))


for fa in range(2, len(sys.argv)):
    fasta = read_fasta(expandvars(sys.argv[fa]))
    ops = editops(target_fasta["fasta"], fasta["fasta"])
    dist = len(ops)

    sample = fasta["name"].split("|")[0].strip()

    org, acc = reference_to_organism_accession[sys.argv[1]]

    results.append(
        {
            "File": os.path.basename(sys.argv[fa]),
            "HasCalls": 0 if dist == 0 else 1,
            "Sample": sample,
            "Organism": org,
            "Accession": acc,
            "Hash": fasta["hash"],
            "Distance": dist,
            # "EditOperations": ops,
        }
    )

with open("/dev/stdout", "w") as f:
    writer = csv.DictWriter(
        f,
        [
            "File",
            "HasCalls",
            "Sample",
            "Organism",
            "Accession",
            "Hash",
            "Distance",
            # "EditOperations"
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(results)
