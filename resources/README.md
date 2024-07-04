This directory is intended for data resources that are imported by the pipeline while it runs. These
data are used as reference materials and consist of a data dump from BOLD (with sequences), a 
release from OpenTree (for topological constraints), and precomputed hidden Markov models for selected
barcode markers. 

The file [bold_ids.tsv](bold_ids.tsv) is a mapping between the taxonomic names connected
to BOLD process IDs and taxa from the Netherlands Species Registry. This file was produced by:

- workflow/scripts/gather_bold_ids.py
- BOLD_Public.05-Apr-2024-curated.tsv
- NSR: https://doi.org/10.15468/rjdpzy (10-06-2024)
