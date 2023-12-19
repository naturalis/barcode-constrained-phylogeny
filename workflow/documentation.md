## Rules

- `create_database` - Puts relevant BOLD data columns into a custom SQLite database.
  Uses BOLD dump TSV as defined in config file. Implemented by the 
  [bold_data_dump.py](scripts/bold_data_dump.py) script.

- `map_opentol` - Enriches the SQLite database with mappings to OpenToL. Because this 
  operates on the same database file as the previous rule, the output is a 0-byte file
  `map_opentol.ok` to indicate that the task was run. Uses the 
  [Open Tree of Life API](https://github.com/OpenTreeOfLife/germinator/wiki/TNRS-API-v3#match_names) 
  to map BOLD taxon names to Open Tree of Life taxonomy IDs. Implemented by the
  [map_opentol.py](scripts/map_opentol.py) script.

### [family_fasta.py](scripts/family_fasta.py)
- Barcodes from the custom database are divided into their taxonomic family groups and written to FASTA files: 'fasta/family/{family}.fasta'

### [hmm_build.py](scripts/hmm_build.py)
- Create a HMM profile for marker.

### [msa_hmm.py](scripts/msa_hmm.py)
- Create MSA using the HMM profile for marker.

### [replace_alignment_ids.py](scripts/replace_alignment_ids.py)
- Replaces the FASTA headers in alignments from '>{barcode_id}' to '>{opentol_id}\_{barcode_id}' 

### [edit_constraint.py](scripts/edit_constraint.py)
- Edits the constraint trees to accomodate for multiple barcodes from one species and prepares the newick files for raxml.

### [distance_matrix.py](scripts/distance_matrix.py)
- Create distance matrix from newick RAXML output.

### [create_submatrices.py](scripts/create_submatrices.py)
- Create submatrix from distance matrix, containing only entries which are found in the opentol database.

### [get_representatives.py](scripts/representatives.py)
- Get the two entries with the highest distance between them. Those are used as representatives.

### [representatives_to_file.py](scripts/representatives_to_file.py)
- Write the two representatives to a fasta file.

### [insert_backbone.py](scripts/insert_backbone.py)
- Insert the subtree into the backbone.


