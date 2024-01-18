# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.0] - 2024-01-11

### Added

- Added the use of a prior for the ambient alpha estimates for RNA and DNA.
  This prevents overfitting and ensures numerical stability at values
  close to 0 or 1. The weight of the prior can be given in the command line
  with the `--amb-prior-w` argument (set to 0.1) by default. The value
  of the prior is updated iteratively in an empirical bayesian fashion,
  in which the average across singlets is taken, weighted by posterior
  probabilities.
- Added an `info` metric for the best assignment in the `summary.txt` file.
  This contains the "effective" number of reads which are weighed by
  their informativeness. Using this, you can approximate the alpha
  value as a beta distribution.

### Changed
- Removed EM iterations for estimation of the droplet ambient fractions
  and fully replaced with maximization of the log marginal likelihood using
  Newton-Raphson. This greatly speeds up model fitting. The Newton-Raphson
  function previously served as the initialization. 
    - Updated the Newton-Raphson function with the prior to ensure that
      the alpha values stay within a numerically stable range.
      There was a bug in some instances when the optimization would occur
      outside of the domain. Fixing this allowed the Newton-Raphson approach
      to converge to the EM estimate, which means the slower EM procedure can be
      removed altogether.

## [0.3.0] - 2023-12-10

### Added

- The summary output file includes `best_singlet` and `best_doublet`
  assignments.

### Changed

- The summary output file directly gives `best_rna_ambient` and
  `best_atac_ambient` estimates. Previously, these had to be calculated
  in a more cubersome manner.
- Improved speed by summing lambda and pi post. probs at the end
  of the Expectation function instead of looping through barcodes again.
- Initialization of alpha is bounded by [0.01, 0.99], changed from [0.05, 0.95].
- More concise verbose output. Include version in output.

### Fixed

- The output for the `best_sample` matches the output of the `best_type`. There
  can be cases where, for example, a singlet sample has the highest likelihood
  among all assignments, but a doublet overall has the highest likelihood
  compared to empty or singlets.
