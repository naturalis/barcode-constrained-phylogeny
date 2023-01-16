import io
import json

import numpy as np
import requests
import os
import argparse
import urllib3
import pandas as pd
import opentree

par_path = os.path.abspath(os.path.join(os.pardir))
http = urllib3.PoolManager()

# User arguments
parser = argparse.ArgumentParser()
parser.add_argument('-outdir1', default=par_path+"/data/",
                    help="Output folder 1: BOLD export data directory")
args = parser.parse_args()

#NOT USED, only fasta
def boldExtract():
    """ Obtains public sequence data from BOLD from the COI-5P marker.

    Downloads records using BOLD's Public Data Portal API. All COI-5P records
     are downloaded and saved as a FASTA file

    Arguments:
        url = String, URL for data retrieval using BOLD's API.
        r: HTTPResponse, variable for retrieving web-url's
        data: response formatted to a string
        fasta_file: file that will hold the data
    """
    # DOES NOT WORK, USE data/bold_COI-5P_Netherlands.fasta
    url = "http://v3.boldsystems.org/index.php/API_Public/sequence?" \
       "marker=COI-5P&taxon=Animalia&geo=Netherlands"
    print('BOLD sequence data retrieval...')
    # Download sequence data from BOLD filtering on marker
    r = requests.get(url)
    # Format byte to string and remove double new lines.
    data = r.content.decode("utf-8").replace('\r\n', '\n')

    # open (and make) fasta file
    fasta_file = open('bold_Netherlands.fasta', 'w')
    # write  string to FASTA file
    fasta_file.write(data)
    # close file
    fasta_file.close()



def boldSpecimen():
    os.makedirs('fasta/family', exist_ok=True)
    url = "http://v3.boldsystems.org/index.php/API_Public/combined?" \
          "&marker=COI-5P&taxon=Animalia&format=tsv"
    print('BOLD specimen and sequence data retrieval...')
    # Download sequence data from BOLD filtering on marker
    r = requests.get(url)
    # Format byte to string and remove double new lines.
    data = r.content.decode('cp1252')
    df = pd.DataFrame([x.split('\t') for x in data.split('\n')[1::]],
                      columns=data.split('\n')[0].split())
    before_filter = len(df)
    df = df[df['marker_codes'].str.contains("COI-5P", na=False)]
    after_filter = len(df)
    print('%i records removed without COI-5P markers. %i barcodes are left.'
          % (before_filter, (before_filter-after_filter)))
    df = df.dropna()
    # loop through unique family names
    for i in set(df['family_name']):
        print('Making FASTA file for sequences from the family %s...' % i)
        # grab records from specific family
        df_order = df[df['family_name'] == i]
        # use function mapOpentol to change species names to Opentol taxonomy
        for species in set(df_order['species_name']):
            y = False
            new_name, y = mapOpentol(species, i, y)
            if y:
                print(df_order.loc[df_order['species_name'] == species])
                df_order['species_name'] = df['species name'].replace([species], new_name)
                print(df_order.loc[df_order['species_name'] == species])



        # make a custom fasta header
        df_order["fasta_header"] = ">" + df_order["processid"].astype(str) + "|" + \
                            df_order["sequenceID"] + "|"+ df_order['species_name']
        # Put header and sequences in fasta format
        fasta_out = df_order['fasta_header'] + "\n" + df_order["nucleotides"]
        pd.set_option('display.expand_frame_repr', False)
        # Save BOLD sequences and their respective header in a fastafile
        np.savetxt("fasta/family/%s.fasta" % i, fasta_out.values,
                   fmt="%s")



def mapOpentol(bold_name, family_name, y):
    import requests
    # Making a POST request
    r = requests.post("https://api.opentreeoflife.org/v3/tnrs/match_names",
                      data=json.dumps({"names": [bold_name],
                                       "context": family_name,
                                       'do_approximate_matching': True,
                                       'include_supressed': False}
                                      )
                      )
    # check status code for response received
    # success code - 200

    if r.ok:
        matches = r.json()['results'][0]['matches']
        if len(matches) != 0:
            matched = matches[0]['matched_name']
            # print content of request
            if matched != bold_name:
                print(matched, 'opentol')
                print(bold_name, 'bold')
                print(matches)
                y = True
                return matched, y
    return bold_name, y



boldSpecimen()
# mapOpentol(['Taschorema nr. apobamum'])