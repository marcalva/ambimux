# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
