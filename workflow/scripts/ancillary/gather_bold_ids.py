import pandas as pd
import requests
import json
import time
from io import StringIO
import sys


def gather_names_NSR(dataframe_taxa):
    """ This gathers the organism names and taxon ID's from a given
    taxonomy file as taken from NSR.
    """
    list_organisms = []
    for index, row in dataframe_taxa.iterrows():
        list_organisms.append(f"{row['genus']} {row['specificEpithet']}")
    list_tax_ids = list(dataframe_taxa['taxonID'])
    return list_organisms, list_tax_ids


def get_bold_organisms(list_organisms):
    """ This looks up organism names in the checklistbank
    to get the organism nomenclature in the BOLD database.
    It returns the organism names in BOLD as a list.
    """
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
        if count % 4000 == 0:
            time.sleep(5)
        count += 1

    return organisms_bold


def get_bold_id(organisms_bold):
    """ This script matches the organism name in the BOLD database to 
    the BOLD IDs that are associated with it.
    It returns a list of BOLD IDs.
    """
    count = 0
    bold_ids = []
    for organism_name in organisms_bold:
        organism_split = organism_name.split(" ")
        try:
            genus = organism_split[0]
            species = organism_split[1]

            response = requests.get(rf"https://v3.boldsystems.org/index.php/API_Public/specimen?taxon={genus}%20{species}&format=tsv")
            res = pd.read_csv(StringIO(response.text), sep="\t")
            bold_ids.append(response.text.split("\n")[1].split("\t")[0])

        except:
            bold_ids.append("")
        count += 1
        if count % 4000 == 0:
            time.sleep(5)
    return bold_ids


if __name__ == '__main__':
    # filename_nsr_taxa = "resources/Taxa.txt"
    filename_nsr_taxa = sys.argv[1]
    dump_file = sys.argv[2] # resources/bold_ids.txt

    dataframe_taxa = pd.read_csv(filename_nsr_taxa)
    list_organisms, list_tax_ids = gather_names_NSR(dataframe_taxa)
    organisms_bold = get_bold_organisms(list_organisms)
    bold_ids = get_bold_id(organisms_bold)

    info_dump = {'NSR_organism': list_organisms[0:len(organisms_bold)],
                    'NSR_ID': list_tax_ids[0:len(organisms_bold)],
                    'BOLD_organism': organisms_bold[0:len(organisms_bold)],
                    'BOLD_ID': bold_ids[0:len(organisms_bold)]} 
    # print(info_dump)
    df_dump = pd.DataFrame(info_dump)
    df_dump.to_csv(dump_file, sep='\t', index=False, header=True)