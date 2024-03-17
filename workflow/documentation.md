## Rules

![](../doc/dag.svg)

---
### `create_database`

Puts relevant BOLD data columns into a custom SQLite database.
Uses BOLD dump TSV as defined in config file. Implemented by the 
[bold_data_dump.py](scripts/create_database.py) script.

---
### `map_opentol` 
 
Enriches the SQLite database with mappings to OpenToL. Because this 
operates on the same database file as the previous rule, the output is a 0-byte file
`map_opentol.ok` to indicate that the task was run. Uses the 
[Open Tree of Life API](https://github.com/OpenTreeOfLife/germinator/wiki/TNRS-API-v3#match_names) 
to map BOLD taxon names to Open Tree of Life taxonomy IDs. Implemented by the
[map_opentol.py](scripts/map_opentol.py) script.

---
### `megatree_loader`

Loads the configured OpenToL topology in the into the SQLite database. This
is a multi-step process that first processes the tree using the conda package
`perl-bio-phylo-forest-dbtree`, exports this from a temporary database to an 
SQL dump, imports that dump into the central database, and then creates a 
0-byte file `megatree_loader.ok`.

---
### `family_fasta`

Barcodes from the SQLite database are divided into their taxonomic family 
groups and written to FASTA files: 'fasta/family/{family}.fasta'. This rule
is run in parallel using `scattergather` and is implemented in the 
[family_fasta.py](scripts/family_fasta.py) script. Note that the parallelization
needs to know ahead of time how many families to expect, which needs to be
specified in the configuration file.

---
### `family_constraint`

For each family FASTA file generates a constraint tree. The constraints are 
generated using the tree indexed by `megatree_loader` and the family FASTA
files from `family_fasta`. This rule is also parallelized using `scattergather`,
and is implemented using the conda package `perl-bio-phylo-forest-dbtree`.

---
### `msa_hmm`

Creates a multiple sequence alignment for each family FASTA file using Hidden
Markov Model alignment. Here, the sequences are also corrected for possible 
revcom issues. It is possible that this is not needed at all because so far 0 
revcom sequences were observed in BOLD. This rule is implemented by the script 
[msa_hmm.py](scripts/msa_hmm.py), which depends on the conda package `hmmer`,
and is parallelized using `scattergather`.

---
### `prep_raxml`

Prepares the aligned ingroup sequences and constraint tree for analysis by 
conda package `raxml-ng`. This step also includes outgroup selection. It is
possible that this rule needs to be reimplemented with an eye on doing 
phylogenetic placement (as opposed to constrained tree searching) instead.
This would obviate the need for outgroup selection and for subsequent re-
rooting. The rule is implemented by the script 
[prep_raxml.py](scripts/prep_raxml.py) and is parallelized using `scattergather`.

---
### `run_raxml`

Runs a phylogenetic analysis for each family FASTA file + constraint tree.
Implemented using the `raxml-ng` conda package and parallelized using
`scattergather`. Might be changed to phylogenetic placement.

---
### `reroot_raxml_output`

Reroots the output from `raxml-ng`. The reason for this is that RAxML emits
unrooted trees when it runs a tree search. The way in which this is done
is by figuring out what the outgroup taxa were and then rooting the tree
on the smallest bipartition that separates the outgroup from the ingroup,
modifying the root branch length by taking the midpoint. There are several
issues with this, because the midpoint is a crude and possibly inappropriate
assumption and because there are cases where selected outgroups nest inside
the ingroup. Implemented by the script 
[reroot_backbone.py](scripts/reroot_backbone.py) and parallelized using
`scattergather`.

---
### `choose_exemplars`

Chooses two exemplars that represent the ingroup in the backbone topology.
At present these are selected by taking the shallowest tip in the ingroup.
This may be inappropriate: perhaps the tallest tips should be selected, or
the median ones (this needs to be tested empirically). Implemented by the
script [choose_exemplars.py](scripts/choose_exemplars.py) and parallelized
using `scattergather`.

---
### `prep_raxml_backbone`

Prepares the backbone constraint tree and alignment by merging the exemplars
from the separate family-level analyses. Implemented by the script
[backbone_constraint.py](scripts/backbone_constraint.py). From this rule
onward, the pipeline is no longer parallel.

---
### `run_raxml_backbone`

Runs the backbone analysis using `raxml-ng`. This is done by a heuristic
search, which may have to be replaced by a phylogenetic placement operation
instead.

---
### `reroot_backbone`

Reroots the backbone tree. This is done by taking the input constraint tree
and reconciling the resulting backbone with the orientation of the input.
Perhaps this is not necessary if we use phylogenetic placement instead.
Implemented by the [reroot_backbone.py](scripts/reroot_backbone.py) script.

---
### `graft_clades`

Grafts the family level subtrees onto the backbone. The branch lengths on
the backbone and those in the family level trees are on a different order
of magnitude, with the (more or less) equivalent branches in the backbone
much longer than in the subtrees. This is presumably because the backbone
spans a broader taxon sample, which consequently has greater saturation and
homoplasy in the alignment, which is accommodated by lengthening branches.
How to scale those relative to one another is a challenge to be addressed
empirically by testing different exemplar selections and rescaling rules.
