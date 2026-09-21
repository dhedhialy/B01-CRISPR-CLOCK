# Exact solver layer

Implementation: `b01_clock/solver/exact.py` + `pruning.py`.

## Object

For ultrametric trees with dated total time \(T\) and global rate bound \([r_{\min}, r_{\max}]\), profile the likelihood over the global rate and compute the profile likelihood-ratio set for each target functional (node age, branch duration, order gap).

## Certificate

Each bound records: topology, rate bounds, anchor \(T\), confidence level, grid resolution, \(L_{\max}\), LR cutoff \(\tfrac12\chi^2_{1,1-\alpha}\), and identification status ∈ {unidentified, partially_identified, point_identified}.

## Validation

Matches finer-grid recomputation within tolerance (`verify_against_grid`). Exhaustive on ≤3 free internal ages.

## Relation to theorem

Layer 1 of the algorithmic architecture: exact enumeration / certified grid optimization on small trees before scalable extensions.
