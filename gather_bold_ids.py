import pandas as pd
import requests
import json
import time
from io import StringIO


def gather_names(dataframe_taxa):
    list_organisms = []
    for index, row in dataframe_taxa.iterrows():
        list_organisms.append(f"{row['genus']} {row['specificEpithet']}")
    return list_organisms


def get_bold_names(list_organisms):
    all_bold_ids = []
    unmatched = []
    count = 0
    length = len(list_organisms)
    write_out_ids = "NSR name\tTaxID\n"
    write_out_umatched = "NSR name\n"
    for organism in list_organisms:
        print(f"{count} out of {length}")
        organism_split = organism.split(" ")
        genus = organism_split[0]
        species = organism_split[1]

        response = requests.get(rf"https://api.checklistbank.org/dataset/37384/nameusage/search?q={genus}%20{species}&offset=0&limit=10")
        try:
            json_resp = response.json()
            bold_id = json_resp['result'][0]['usage']['name']['link']
            # "http://bins.boldsystems.org/index.php/Taxbrowser_Taxonpage?taxid=12439"
            bold_id = bold_id.split("taxid=")[1]
            all_bold_ids.append(bold_id)
            write_out_ids += f"{genus} {species}\t{bold_id}\n"
        except:
            unmatched.append(organism)
            write_out_ids += f"{genus} {species}\n"
        if count % 4000 == 0:
            time.sleep(5)
        count += 1
    
    open("resources/bold_ids.txt", "w").write(write_out_ids)
    open("resources/unmatched.txt", "w").write(write_out_umatched)


def match_unmatched(filename="resources/unmatched.txt"):
    list_organisms =  open(filename, "r").readlines()[1:]
    all_bold_ids = []
    unmatched = []
    count = 0
    length = len(list_organisms)
    write_out_ids = ""
    write_out_umatched = "NSR name\n"
    for organism in list_organisms:
        print(f"{count} out of {length}")
        organism_split = organism.split(" ")
        genus = organism_split[0]
        species = organism_split[1]

        response = requests.get(rf"https://v3.boldsystems.org/index.php/API_Public/specimen?taxon={genus}%20{species}&format=tsv")
        try:
            gathered_ids = set(pd.read_csv(StringIO(response.text), sep="\t")["species_taxID"].tolist())
            prepped_id = str(int(list({x for x in gathered_ids if pd.notna(x)})[0]))
            all_bold_ids.append(prepped_id)
            write_out_ids += f"{genus} {species}\t{prepped_id}\n"
        except:
            unmatched.append(organism)
            write_out_umatched += f"{genus} {species}\n"
        if count % 4000 == 0:
            time.sleep(5)
        count += 1
    open("resources/bold_names.txt", "a").write(write_out_ids)
    open("resources/unmatched.txt", "w").write(write_out_umatched)


if __name__ == '__main__':
    filename_nsr_taxa = "resources/Taxa.txt"
    dataframe_taxa = pd.read_csv(filename_nsr_taxa)
    list_organisms = gather_names(dataframe_taxa)
    get_bold_names(list_organisms)
    match_unmatched()
