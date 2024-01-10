Bactria pipeline results
========================

This folder structure contains the output of the pipeline in accordance with the
configuration specified in the corresponding code release, which in a local install
can be found in ../config/config.yaml

The subfolder structure is organized as follows:

- [`databases`](databases) contains the aggregated SQLite database with barcodes,
  BOLD taxonomy information reconciled with the Open Tree of Life, and the topology
  of the open tree. Any zero-byte files with the \*.ok extension that are in this 
  folder can be ignored as they functioned as status updates for the various 
  database update processes.
- [`fasta`](fasta) principally contains FASTA multiple sequence files (\*.fa) and
  Newick tree files (\*.tre). These are inputs for RAxML, which produces a number
  of output files, all of which will have \*.raxml.\* in their file names. The files
  are organized in families, which are processed in parallel and which are organized
  in folders with m-of-n names. The backbone files across families are in the root
  of the `fasta` folder.

