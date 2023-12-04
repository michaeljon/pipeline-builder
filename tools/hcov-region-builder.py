#!/usr/bin/env python

import sys
import json

feature_data = {}
organisms = {
    # HCoV
    "AF304460.1": "hcov-229e",
    "AY597011.2": "hcov-hku1",
    "AY567487.2": "hcov-nl63",
    "AY585228.1": "hcov-oc43",
    # SARS-CoV-2
    "MN908947.3": "sars-cov-2",
    # HRV
    "NC_001617.1": "hrv-a",
    "NC_038311.1": "hrv-a1",
    "NC_001490.1": "hrv-b14",
    "NC_038312.1": "hrv-b3",
    "NC_009996.1": "hrv-c",
    "NC_038878.1": "hrv-nat001",
}

for i, (k, v) in enumerate(organisms.items()):
    feature_data[v] = {"sequence": k, "name": v, "region": {"start": -1, "stop": -1}, "genes": {}}

for i, f in enumerate(sys.argv):
    if i > 0:
        with open(f, "r") as gff:
            while line := gff.readline().rstrip():
                if line.startswith("#") == False:
                    elements = line.split("\t")
                    sequence_id = organisms[elements[0]]
                    start_pos = int(elements[3])
                    stop_pos = int(elements[4])

                    if elements[2] == "region":
                        feature_data[sequence_id]["region"]["start"] = start_pos
                        feature_data[sequence_id]["region"]["stop"] = stop_pos

                    if elements[2] == "five_prime_UTR":
                        feature_data[sequence_id]["genes"]["5-prime-utr"] = {"start": start_pos, "stop": stop_pos}

                    if elements[2] == "three_prime_UTR":
                        feature_data[sequence_id]["genes"]["3-prime-utr"] = {"start": start_pos, "stop": stop_pos}

                    elif elements[2] == "gene":
                        infos = elements[8].split(";")
                        if infos != None and len(infos) > 0:
                            for i in infos:
                                v = i.split("=")
                                if i.startswith("Name="):
                                    name = v[1]

                                    feature_data[sequence_id]["genes"][name] = {"start": start_pos, "stop": stop_pos}


print(json.dumps(feature_data, indent=2, sort_keys=False))
