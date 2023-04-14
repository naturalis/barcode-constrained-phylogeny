# barcode-constrained-phylogeny
This repository contains code and data for building topologically-constrained phylogenies.

The international database BOLD Systems contains DNA barcodes for hundreds of thousands of species, with multiple barcodes per species. Theoretically, this data could be filtered and aligned per DNA marker to make phylogenetic trees. However, there are two limiting factors: there is a maximum number and size of trees that can be built, and barcodes are not considered ideal for building large trees because they do not provide a strong signal.

Both problems can be solved by using the Open Tree of Life as a backbone. The BOLD data can be split into chunks that correspond to the Open Tree of Life. These chunks can be made into alignments and trees, and the backbone can also be used as a constraint in the algorithms to create these.

This repository will prototype the use of Open Tree of Life as a backbone for both animal and plant species.

## Scripts
### [unzip_targz.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/unzip_targz.py)
- Unzips a targz file, more specifically a [datarelease](https://www.boldsystems.org/index.php/datapackage?id=BOLD_Public.30-Dec-2022) from BOLD Systems containing a snapshot of the barcode database (more than 8 million barcodes as of 30-DEC-2022).

### [bold_data_dump.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/bold_data_dump.py) 
- Puts relevant BOLD data columns into a custom SQLite database.

### [alter_tables.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/alter_tables.py)
- Manipulates the BOLD data in the custom database and makes two tables, one for taxon data and one for barcode entries.

### [map_opentol.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/map_opentol.py)
- Uses the [Checklisbank API](https://api.checklistbank.org/) to map BOLD taxon names to Open Tree of Life taxonomy IDs. 

### [family_fasta.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/family_fasta.py)
- Barcodes from the custom database are divided into their taxonomic family groups and written to FASTA files: 'fasta/family/{familyname}.fasta'

### [download_macse.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/macse/src/download_macse.py)
- Download program macse (Only works on linux distributions).

### [create_MSA](https://github.com/naturalis/barcode-constrained-phylogeny/blob/macse/src/create_MSA.py)
- Using macse create a Multiple Sequence Alignment in files: '{filename_NT}.fasta' and '{filename_AA}.fast'.

### [create_tree](https://github.com/naturalis/barcode-constrained-phylogeny/blob/macse/src/create_tree.py)
- A distance matrix is made based on a NT/AA file. From the distance matrix a UPGMA tree is made. 
