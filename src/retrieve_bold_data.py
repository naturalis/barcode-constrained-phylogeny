import json
import numpy as np
import os
import argparse
import urllib3
import pandas as pd
import requests

# Directory
par_path = os.path.abspath(os.path.join(os.pardir))
http = urllib3.PoolManager()

# User arguments
parser = argparse.ArgumentParser()
parser.add_argument('-outdir1', default=par_path+"/data/",
                    help="Output folder 1: BOLD export data directory")
args = parser.parse_args()


def boldSpecimen():
    """ Obtains public barcode sequence data from BOLD from the COI-5P marker.
    Download records using BOLD's Public Data Portal API. All COI-5P records
    are downloaded. The species names from BOLD are searched for OpenTol
    taxonomy using the function maptoOpentol(). When another species name is
    used in OpenTol it gets replaced in the DataFrame containing the BOLD
    records. A fasta file per family group is made using custom metadata and
    saved in the folder 'fasta/family/{family_name}'

    Arguments:
        url = String, URL for data retrieval using BOLD's API.
        r: HTTPResponse, variable for retrieving web-url's
        data: response formatted to a string
        df: data frame containing DNA barcode records from BOLD.
        df_family: dataframe with the BOLD records from one family group.
        fasta_out: BOLD data reformatted to FASTA format
    """
    # Make directory if it does not exist yet.
    os.makedirs('fasta/family', exist_ok=True)
    # BOLD url for specimen records with arguments
    url = "http://v3.boldsystems.org/index.php/API_Public/combined?" \
          "&marker=COI-5P&taxon=Animalia&format=tsv"
    print('BOLD specimen and sequence data retrieval...')
    # Download sequence data from BOLD filtering on marker and taxon
    r = requests.get(url)
    # Format byte to string and remove double new lines.
    data = r.content.decode('cp1252')
    df = pd.DataFrame([x.split('\t') for x in data.split('\n')[1::]],
                      columns=data.split('\n')[0].split())
    df = df.dropna()
    pd.set_option('display.max_columns', None)
    print(df.tail(10))
    before_filter = len(df)
    # Filter records on COI-5P marker
    df = df[df['marker_codes'].str.contains("COI-5P", na=False)]
    after_filter = len(df)
    print(before_filter , "before")
    print(after_filter , "after COI-5P filter")

    id = 0
    # Loop through unique family (taxonomy) names
    for family in set(df['family_name']):
        print('Making FASTA file for sequences from the family %s...' % family)
        # Grab records from specific family
        df_family = df[df['family_name'] == family]

        # Use function mapOpentol to change species names to Opentol taxonomy
        for species_name in set(df_family['species_name']):
            id += 1
            opentol_name, found_new = mapOpentol(species_name, family)
            if found_new:
        # Disable panda warnings
                pd.set_option('display.expand_frame_repr', False)
                pd.options.mode.chained_assignment = None

        # Make a custom fasta header
                sub_df = df_family[df_family['species_name'] == species_name]
                sub_df["fasta_header"] = ">" + df['species_name'].replace(
                    [species_name], (str(opentol_name) + "|" + str(id + 1)))
                print(sub_df)
    # +\
    # "|" + df_family["sequenceID"] + \
        # Put header and sequences in fasta format
                fasta_out = sub_df['fasta_header'] + "\n" + sub_df["nucleotides"]

                # Save BOLD sequences and their respective header in a fastafile
                np.savetxt("fasta/family/%s.fasta" % family, fasta_out.values,
                           fmt="%s")



def mapOpentol(species_name, family):
    """ Searches OpenTol taxonomy for a match with the BOLD species name.
    Search constraints are BOLD species name, family name the species belongs
    to, fuzzy matching enabled and filtered on flagged on not accepted species.
        Arguments:
            species_name: String containing species name from BOLD.
            family: String containing corresponding family name from BOLD.
            url = String, URL for data retrieval using OpenTOL's API.
            r: HTTPResponse, variable for retrieving web-url's.
            matches: response formatted to a list containing the matches.
            match: String with matched OpenTOL species name.
        Returns:
            species_name: String containing old or new species name.
            boolean: True if the name is updated, False if not
        """
    # Making a POST request with fuzzy matching(approximate matching)
    r = requests.post("https://api.opentreeoflife.org/v3/tnrs/match_names",
                      data=json.dumps({"names": [species_name],
                                       "context": family,
                                       'do_approximate_matching': True,
                                       'include_supressed': False}
                                      )
                      )
    # Check if http response is valid
    if r.ok:
        # Grab the matched species names from json result
        matches = r.json()['results'][0]['matches']
        # Check if there are any matches
        if len(matches) != 0:
            # Grab the OpenTol species name that matched to the BOLD one.
            matched = matches[0]['taxon']['ott_id']
            # Check if they are the same
                # If not return the new name and the boolean True
            return matched, True
    # Return the original species name and False if there is no match
    return species_name, False


if __name__ == '__main__':
    boldSpecimen()
