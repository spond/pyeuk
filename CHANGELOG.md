# Changelog

## 0.4.0

Minor rather than patch: the default distance weighting changed and the public API grew.
Callers who relied on `0.3.0`'s weighting will get different distances from the same sheet.

### Changed

- **`PyEukDistanceEngine` default weighting is now `weight_mode="heterozygosity"`** (`2p(1-p)`)
  rather than the KING form `1/sqrt(p(1-p))` (#11). The KING form is a standardisation for
  centred genotype dosages; applied to presence/absence indicator columns it is inverted
  relative to intent, being minimised at `p = 0.5` and diverging as `p -> 0`, so it weights
  rare alleles hardest while the discriminating signal sits in balanced columns. Pass
  `weight_mode="king"` to restore the previous behaviour.
- `CyclosporaClusterFinder.find_clusters` now honours `k_min` in label-free mode (#12).
  Previously `k_min` was applied only in the supervised branch and silently ignored in the
  default one, so a caller who knew the cohort contained three groups had no way to say so.
- The `k = 1` fallback message now reports which criterion each candidate failed. It
  previously named the minimum-cluster-size guard unconditionally, including when the
  relative-gap floor was what rejected the candidate.

### Added

- `weight_mode`, `min_maf` and `project_psd` arguments on the distance engine (#11).
- `relative_gap_floor` as a `CyclosporaClusterFinder` constructor argument and a
  `find_clusters` per-call override. The `0.2200` threshold was previously a hardcoded
  literal, which is a reasonable default and an awkward only-value: the appropriate floor
  depends on how many groups are expected, since the more real groups a cohort has, the
  smaller each merge-height gap is relative to the root.
- `CyclosporaClusterFinder.last_selection_meta`, a dict describing how `k` was chosen.
  `status` is one of `optimal`, `floor_override`, `unsatisfiable_constraint` or
  `single_group`. The last two both return `k = 1`, and the field is what lets a caller
  distinguish a genuinely single-group cohort from an unsatisfiable `k_min`.

### Documentation

- `PROVENANCE.md` for the Cryptosporidium and Giardia example data now states that the
  fixtures are synthetic, and the generating script and README agree with it (#10).
