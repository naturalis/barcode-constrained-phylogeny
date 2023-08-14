## Scripts
### [unzip_targz.py](unzip_targz.py)
- Unzips a targz file, more specifically a [datarelease](https://www.boldsystems.org/index.php/datapackage?id=BOLD_Public.30-Dec-2022) from BOLD Systems containing a snapshot of the barcode database (more than 8 million barcodes as of 30-DEC-2022).

### [bold_data_dump.py](bold_data_dump.py) 
- Puts relevant BOLD data columns into a custom SQLite database.

### [map_opentol.py](map_opentol.py)
- Uses the [Open Tree of Life API](https://github.com/OpenTreeOfLife/germinator/wiki/TNRS-API-v3#match_names) to map BOLD taxon names to Open Tree of Life taxonomy IDs. 

### [family_fasta.py](family_fasta.py)
- Barcodes from the custom database are divided into their taxonomic family groups and written to FASTA files: 'fasta/family/{family}.fasta'

### [create_MSA.py](create_MSA.py)
- Uses masce to create a Multiple Sequence Alignment as output files:'{family_NT}.fasta' and '{family_AA}.fasta'.

### [replace_alignment_ids.py](replace_alignment_ids.py)
- Replaces the FASTA headers in alignments from '>{barcode_id}' to '>{opentol_id}\_{barcode_id}' 

### [edit_constraint.py](edit_constraint.py)
- Edits the constraint trees to accomodate for multiple barcodes from one species and prepares the newick files for raxml.

### [distance_matrix.py](distance_matrix.py)
- Create distance matrix from newick RAXML output.

### [create_submatrices.py](create_submatrices.py)
- Create submatrix from distance matrix, containing only entries which are found in the opentol database.

### [get_representatives.py](representatives.py)
- Get the two entries with the highest distance between them. Those are used as representatives.

### [representatives_to_file.py](representatives_to_file.py)
- Write the two representatives to a fasta file.

### [insert_backbone.py](insert_backbone.py)
- Insert the subtree into the backbone.


