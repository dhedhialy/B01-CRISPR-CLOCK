# ConvexML: Scalable and accurate inference of single-cell chronograms from CRISPR/Cas9 lineage tracing data

[https://doi.org/10.5061/dryad.qrfj6q5nz](https://doi.org/10.5061/dryad.qrfj6q5nz)

Our simulated trees paired with lineage tracing data encompass a large number of lineage tracing regimes, which are used to assess the performance of our proposed branch length estimator.

## Description of the data and file structure

For each lineage tracing regime, 50 simulations are performed. All trees have exactly 400 leaves, and were simulated as described in the manuscript. The `default' regime consists of:

13 barcodes.

3 target sites per barcode.

mutation rate adjusted to obtain an expected 50% mutated entries in the character matrix.

100 indel states.

20% missing data, with 10% coming from heritable epigenetic silencing and 10% coming from sequencing dropouts. (This does not include missing data further introduced by double-resection events, which we also simulate.)

Each lineage tracing regime is obtained by perturbing this 'default' lineage tracing regime by varying one of the above parameters. Specifically, we consider varying:

* number_of_cassettes: number of barcodes (a.k.a. cassettes) in the set {3, 6, 13, 20, 30, 50} (with 13 being the default)
* number_of_states: number of states in the set {5, 10, 25, 50, 100, 500, 1000} (with 100 being the default)
* expected_proportion_mutated: expected proportion mutated in the set {10%, 30%, 50%, 70%, 90%} (with 50 being the default)
* expected_prop_missing: percent missing from epigenetic silencing and sequencing dropouts in the set {0%, 10%, 20%, 30%, 40%, 50%, 60%}, with the percent coming from sequencing dropouts fixed to 10% (except when the total is 0%, in which case it is set to 0%)

The data from each simulation is stored specifying the parameter that was varied, so, for example, the simulated data when the number of barcodes is 30 is stored under "trees/number_of_cassettes/30/". In this directory, for each repetition, we have three files:

* tree_{repetition}_character_matrix.csv: Contains the lineage tracing data in CSV format.
* tree_{repetition}_newick.txt: Contains the tree in newick format, with branch lengths.
* tree_{repetition}_fitness.txt: Contains the fitness of each leaf node, with one line per leaf, containing the leaf id and its fitness.
* tree_{repetition}_CassiopeiaTree.pkl: Contains the pickled CassiopeiaTree object from the simulation, which in particular contains the fitness of different nodes in the tree, ancestral lineage tracing barcodes, etc. It is not necessary for reproducing any of our results, but we provided it in case it is convenient..

## Code/Software

Trees and lineage tracing data were simulated using the Cassiopeia open source Python package, as described in our manuscript.