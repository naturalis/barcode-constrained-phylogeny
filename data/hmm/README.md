# Hidden Markov models

This folder contains Hidden Markov Models as produced by `hmmbuild`. The purpose of 
these is twofold:

1. As a quick way to check the orientation of sequences, because sometimes (possibly)
   we encounter reverse-complemented (3') sequences when we always want the other
   orientation. This can be checked with `hmmalign`, because the alignment with the
   greater number of maximum posterior probability columns (indicated by `*` in a 
   Stockholm alignment) must be the right orientation.
2. To align the sequences linearly, i.e. without the need for multiple sequence
   alignment of chunks. We can simply align each individually in O(n) time.

These files follow exact naming conventions. Everything before the `.hmm` must be the
name of a marker as defined in configuration file.