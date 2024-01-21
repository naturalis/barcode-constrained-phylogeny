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

## Files contained and planned

- The file [COI-5P.hmm](COI-5P.hmm) was obtained from the [FinPROTAX](https://github.com/psomervuo/FinPROTAX) project,
  and corresponds to the file `FinPROTAX/modelCOIfull/refs.hmm` contained within the ZIP file made available from their
  project repository.
- [matK_rbcL.hmm](matK_rbcL.hmm) was produced by Noah Scheffer from a concatenated alignment of matK and rbcL produced
  during her internship. *TODO* locate that alignment and add it here for reference.
- For ITS, it seems like it's possible to create an HMM by taking the following steps:
  1. Create a reference alignment from taxonomically broadly sampled sequences. This alignment should probably be done
     using `mafft` and should be carefully curated and checked by eye. The alignments tend to be extremely gappy but 
     also have conserved regions, which we will focus on. Those regions should be identifiable by their relative 
     coverage in the alignment and should be extracted, e.g. by a rule that retains all columns with >50% coverage.
  2. Build an HMM from the pruned alignment with the conserved regions using `hmmbuild`.
  3. Run sequences against the alignment using `hmmalign --trim`. These too will be very gappy. Hopefully it is possible
     to figure out which of the columns correspond with those in the pruned HMM, and which ones were inserted. Prune
     the ones that were inserted.