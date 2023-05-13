# barcode-constrained-phylogeny
This repository contains code and data for building very large, topologically-constrained barcode phylogenies through a divide-and-conquer strategy. Such trees are useful as reference materials in the comparable calculation of alpha and beta biodiversity metrics across metabarcoding assays. The input data for the approach we develop here comes from BOLD. The international database [BOLD Systems](https://www.boldsystems.org/index.php) contains DNA barcodes for hundred of thousands of species, with multiple barcodes per species. Theoretically, this data could be filtered and aligned per DNA marker to make phylogenetic trees. However, there are two limiting factors: building very large phylogenies is computationally intensive, and barcodes are not considered ideal for building big trees because they are short (providing insufficient signal to resolve large trees) and because they tend to saturate across large patristic distances.

Both problems can be mitigated by using the [Open Tree of Life](https://tree.opentreeoflife.org/opentree/argus/opentree13.4@ott93302) as a further source of phylogenetic signal. The BOLD data can be split into chunks that correspond to Open Tree of Life clades. These chunks can be made into alignments and subtrees. The OpenTOL can be used as a constraint in the algorithms to make these. The chunks are then combined in a large synthesis by grafting them on a backbone made from exemplar taxa from the subtrees. Here too, the OpenTOL is a source of phylogenetic constraint.

In this repository this concept is prototyped for both animal species and plant species.

## Installation

The pipeline and its dependencies are managed using conda. On a linux or osx system, you can follow these steps to set up the `bactria` Conda environment using an `environment.yml` file and a `requirements.txt` file:

1. **Clone the Repository:**  
   Clone the repository containing the environment files to your local machine:
   ```bash
   git clone https://github.com/naturalis/barcode-constrained-phylogeny.git
   cd barcode-constrained-phylogeny
2. **Create the Conda Environment:**
   Create the bactria Conda environment using the environment.yml file with the following command:
   ```bash
   conda env create -f environment.yml
   ```
   This command will create a new Conda environment named bactria with the packages specified in the environment.yml file. This file also includes pip packages specified in the requirements.txt file, which will be installed after the Conda packages.
3. **Activate the Environment:**
   After creating the environment, activate it using the conda activate command:
   ```bash
   conda activate bactria
   ```
4. **Verify the Environment:**
   Verify that the bactria environment was set up correctly and that all packages were installed using the conda list command:
   ```bash
   conda list
   ```
   This command will list all packages installed in the active Conda environment. You should see all of the packages specified in the environment.yml file and the requirements.txt file.

## How to run

**Explain here how to run the snakemake targets**

## Repository layout

All data used and generated are located in the [data/](https://github.com/naturalis/barcode-constrained-phylogeny/tree/main/data/) directory. 
The snakefile and python scripts are found in the [src/](https://github.com/naturalis/barcode-constrained-phylogeny/tree/main/src/) directory. 


## Scripts
### [unzip_targz.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/unzip_targz.py)
- Unzips a targz file, more specifically a [datarelease](https://www.boldsystems.org/index.php/datapackage?id=BOLD_Public.30-Dec-2022) from BOLD Systems containing a snapshot of the barcode database (more than 8 million barcodes as of 30-DEC-2022).

### [bold_data_dump.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/bold_data_dump.py) 
- Puts relevant BOLD data columns into a custom SQLite database.

### [map_opentol.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/map_opentol.py)
- Uses the [Open Tree of Life API](https://github.com/OpenTreeOfLife/germinator/wiki/TNRS-API-v3#match_names) to map BOLD taxon names to Open Tree of Life taxonomy IDs. 

### [family_fasta.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/family_fasta.py)
- Barcodes from the custom database are divided into their taxonomic family groups and written to FASTA files: 'fasta/family/{family}.fasta'

### [create_MSA](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/create_MSA.py)
- Uses masce to create a Multiple Sequence Alignment as output files:'{family_NT}.fasta' and '{family_AA}.fasta'.

### [replace alignment ids](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/replace_alignment_ids.py)
- Replaces the FASTA headers in alignments from '>{barcode_id}' to '>{opentol_id}\_{barcode_id}' 

### [edit_constraint](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/edit_constraint.py)
- Edits the constraint trees to accomodate for multiple barcodes from one species and prepares the newick files for raxml.


## License

MIT License

Copyright (c) 2022 Naturalis Biodiversity Center

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
