# barcode-constrained-phylogeny
This repository contains code and data for building topologically-constrained phylogenies.

The internationsl database [BOLD Systems](https://www.boldsystems.org/index.php) contains DNA barcodes for hundred of thousands of species, with multiple barcodes per species. Theoretically, this data could be filtered and aligned per DNA marker to make phylogenetic trees with. However, there are two limiting factors: there is a maximum of the number and size trees that can be build, and barcodes are not considered ideal for building big trees because they do not give a strong signal. 

Both problems can be solved by using the [Open Tree of Life](https://tree.opentreeoflife.org/opentree/argus/opentree13.4@ott93302) as a backbone. The BOLD data can be split into chunks that are corresponding to Open Tree of Life. These chunks can be made into alignments and trees. The backbone can also be used as a constraint in the algorithms to make these.

In this repository this will be prototyped for both animal species and plant species.

## Scripts
### [unzip_targz.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/unzip_targz.py)
- Unzips a targz file, more specifically a [datarelease](https://www.boldsystems.org/index.php/datapackage?id=BOLD_Public.30-Dec-2022) from BOLD Systems containing a snapshot of the barcode database (more than 8 million barcodes as of 30-DEC-2022).

### [bold_data_dump.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/bold_data_dump.py) 
- Puts relevant BOLD data columns into a custom SQLite database.


### [alter_tables.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/alter_tables.py)
- Manipulates the BOLD data in the custom database and makes two tables, one for taxon data and one for barcode entries.


### [map_to_opentol.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/map_opentol.py)
- Uses the [Checklisbank API](https://api.checklistbank.org/) to map BOLD taxon names to Open Tree of Life taxonomy IDs. 
