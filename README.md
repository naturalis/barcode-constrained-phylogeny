# barcode-constrained-phylogeny
This repository contains code and data for building very large, topologically-constrained barcode phylogenies through a divide-and-conquer strategy. Such trees are useful as reference materials in the comparable calculation of alpha and beta biodiversity metrics across metabarcoding assays. The input data for the approach we develop here comes from BOLD. The international database [BOLD Systems](https://www.boldsystems.org/index.php) contains DNA barcodes for hundred of thousands of species, with multiple barcodes per species. Theoretically, this data could be filtered and aligned per DNA marker to make phylogenetic trees. However, there are two limiting factors: building very large phylogenies is computationally intensive, and barcodes are not considered ideal for building big trees because they are short (providing insufficient signal to resolve large trees) and because they tend to saturate across large patristic distances.

Both problems can be mitigated by using the [Open Tree of Life](https://tree.opentreeoflife.org/opentree/argus/opentree13.4@ott93302) as a further source of phylogenetic signal. The BOLD data can be split into chunks that correspond to Open Tree of Life clades. These chunks can be made into alignments and subtrees. The OpenTOL can be used as a constraint in the algorithms to make these. The chunks are then combined in a large synthesis by grafting them on a backbone made from exemplar taxa from the subtrees. Here too, the OpenTOL is a source of phylogenetic constraint.

In this repository this concept is prototyped for both animal species and plant species.

## Requirements/dependencies

- [raxml-ng v1.1.0](https://github.com/amkozlov/raxml-ng/releases)
- [macse v2.06](https://bioweb.supagro.inra.fr/macse/releases/macse_v2.06.jar)
- [sqlite3](https://sqlite.org/download.html) **which version?**
- python **which version?**

Further dependencies are specified in [requirements.txt](requirements.txt)

## How to install

FILL ME

## How to run

FILL ME

## Repository layout

FILL ME

## Scripts
### [unzip_targz.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/unzip_targz.py)
- Unzips a targz file, more specifically a [datarelease](https://www.boldsystems.org/index.php/datapackage?id=BOLD_Public.30-Dec-2022) from BOLD Systems containing a snapshot of the barcode database (more than 8 million barcodes as of 30-DEC-2022).

### [bold_data_dump.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/bold_data_dump.py) 
- Puts relevant BOLD data columns into a custom SQLite database.

### [alter_tables.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/alter_tables.py)
- Manipulates the BOLD data in the custom database and makes two tables, one for taxon data and one for barcode entries.

### [map_opentol.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/map_opentol.py)
- Uses the [Open Tree of Life API](https://github.com/OpenTreeOfLife/germinator/wiki/TNRS-API-v3#match_names) to map BOLD taxon names to Open Tree of Life taxonomy IDs. 

### [family_fasta.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/family_fasta.py)
- Barcodes from the custom database are divided into their taxonomic family groups and written to FASTA files: 'fasta/family/{familyname}.fasta'

### [download_macse.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/macse/src/download_macse.py)
- Download program macse (Only works on linux distributions).

### [create_MSA](https://github.com/naturalis/barcode-constrained-phylogeny/blob/macse/src/create_MSA.py)
- Using macse create a Multiple Sequence Alignment in files: '{filename_NT}.fasta' and '{filename_AA}.fast'.

### [create_tree](https://github.com/naturalis/barcode-constrained-phylogeny/blob/macse/src/create_tree.py)
- A distance matrix is made based on a NT/AA file. From the distance matrix a UPGMA tree is made. 

## License

FILL ME
