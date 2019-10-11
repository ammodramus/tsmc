## Triploid SMC (TSMC) demographic inference method

This repository contains a one-off inference method for inferring the evolutionary history of asexual triploid individuals who have most recently undergone [apomictic](https://en.wikipedia.org/wiki/Apomixis) reproduction for some number of generations following an ancient transition from diploid sexual reproduction to triploid apomictic reproduction. The method takes as input a sequence of 0's, 1's and 2's representing the three possible (unphased) triploid genotypes in each bin of some number of base-pairs, similar to the original PSMC method. It infers the time since the transition from diploid sexuality to triploid asexuality together with the population-size history of the diploid ancestors.

The model underlying the inference method is a sequentially Markov coalescent model of three unphased chromosomes, where the state at each position in the genome is the vector (t_3, t_2) containing the more recent and more ancient coalescence times between the three chromosomes of the triploid genome. Because we model unphased genotypes rather than phased haplotypes, it is unnecessary to model the topology of the tree at each position.

To compile, type `make` from the `src/` subdirectory. Requires gcc 4+. See `tsmc -h` for usage and [`src/tsmcequations.pdf`] for model details.

The method was originally intended to be used as part of the *Potamopyrgus antipodarum* genome project, as this species features many independent transitions from diploid sexual ancestors to polyploid apomictic descdendant lineages. However, the *P. antipodarum* genome has also undergone a very recent duplication, which has thus far prevented straightforward application of this method.
