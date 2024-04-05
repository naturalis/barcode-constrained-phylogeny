This directory develops a test to see which exemplar picking strategy is best. The analysis
is implemented in a separate conda environment, i.e.:

```{bash}
mamba env create -f env.yml
mamba activate exemplars_test
```

Then, to set things up, the shell script [estimate.sh](estimate.sh) is run. This script
iterates over the three trees in this folder. The trees result from the pipeline being
run for the Primates, using the three picking strategies (i.e. picking the tallest, shortest,
or median tips in the clades on either side of the root of the subtrees). The script fetches
the corresponding sequences from the databases, aligns them with `hmmalign` and then estimates
the branch lengths with `raxml`. This procedure is a bit sketchy in that, unlike in the 
pipeline, here we don't check for reverse complemented sequences. This means that a small 
number of sequences will be horrendously aligned (about 3-4 out of 300+). We'll remove these
later. It doesn't matter.

The script is just run 'as is', i.e.:

```{bash}
sh estimate.sh
```

Note that this will produce a bunch of files inside this directory, as `raxml` is wont to do.

Then, the next steps are done by R, with the [exemplars.Rmd](exemplars.Rmd) notebook.
