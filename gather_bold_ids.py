import pandas as pd
import requests
import json
import time
from io import StringIO

def gather_names_NSR(dataframe_taxa):
    list_organisms = []
    for index, row in dataframe_taxa.iterrows():
        list_organisms.append(f"{row['genus']} {row['specificEpithet']}")
    list_tax_ids = list(dataframe_taxa['taxonID'])
    return list_organisms, list_tax_ids


def get_bold_organisms(list_organisms):
    count = 0
    length = len(list_organisms)
    organisms_bold = []
    for organism in list_organisms:
        # print(f"{count} out of {length}")
        organism_split = organism.split(" ")
        genus = organism_split[0]
        species = organism_split[1]

        # checklist bank api
        response = requests.get(rf"https://api.checklistbank.org/dataset/37384/nameusage/search?q={genus}%20{species}&offset=0&limit=10")
        try:
            json_resp = response.json()
            organism_name = json_resp['result'][0]['usage']['name']['scientificName']
            # "http://bins.boldsystems.org/index.php/Taxbrowser_Taxonpage?taxid=12439"
            organisms_bold.append(organism_name)
        except:
            organisms_bold.append(f"{genus} {species}")
        count += 1

    return organisms_bold


def match_unmatched(organisms_bold):
    # lines = open(filename, "r").readlines()[1:]
    count = 0
    # length = len(lines)
    bold_ids = []
    for organism_name in organisms_bold:
        # print(f"{count} out of {length}")
        organism_split = organism_name.split(" ")
        try:
            genus = organism_split[0]
            species = organism_split[1]

            response = requests.get(rf"https://v3.boldsystems.org/index.php/API_Public/specimen?taxon={genus}%20{species}&format=tsv")
            res = pd.read_csv(StringIO(response.text), sep="\t")
            gathered_ids = set(res["processid"].tolist())
            prepped_id = str(list({x for x in gathered_ids if pd.notna(x)})[0])
            bold_ids.append(prepped_id)
        except:
            bold_ids.append("")
        count += 1
    return bold_ids


if __name__ == '__main__':
    filename_nsr_taxa = "resources/Taxa.txt"
    dataframe_taxa = pd.read_csv(filename_nsr_taxa)
    list_organisms, list_tax_ids = gather_names_NSR(dataframe_taxa)
    organisms_bold = get_bold_organisms(list_organisms)
    bold_ids = match_unmatched(organisms_bold)

    info_dump = {'NSR_organism': list_organisms[0:len(organisms_bold)],
                    'NSR_ID': list_tax_ids[0:len(organisms_bold)],
                    'BOLD_organism': organisms_bold[0:len(organisms_bold)],
                    'BOLD_ID': bold_ids[0:len(organisms_bold)]} 
    # print(info_dump)
    df_dump = pd.DataFrame(info_dump)
    df_dump.to_csv("resources/bold_ids.txt", sep='\t', index=False, header=True)
