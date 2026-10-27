# Roadmap

Planned work after the restructuring of v2 (pull request #6). Version 2 has the same results
as the production pipeline; the iterations below may change that, each in its own branch and
tested on Nibbler against the production baseline (see `docs/nibbler-test.md`).

Open decisions about current behaviour stay in `docs/known-issues.md`; this file lists work to
do. Remove an item when it is merged.

## Small follow-ups

- **D14: `isec_to_vcf.py` fails on contigs other than chr1-22, X and Y.** Sort other contigs
  after `chrY`, by name. See `docs/known-issues.md`.
- **Make `CircuSNV.sh` executable in git** (`git update-index --chmod=+x CircuSNV.sh`). It is
  stored as 100644, so a `chmod +x` on the cluster shows up as a local change and blocks
  `git pull`.
- **Extend `.gitignore`** with the files that live in the installation directory but are not
  in git: `env/`, `MNVista`, `MNV_PoN_No_overlap_new.txt`, `hg38_genome.UNMET_bin10bp.bed`,
  `PoN/`, `CircuSNV.out`, `CircuSNV.err`, and test output directories (`new_output*/`,
  `old_output*/`). The repository is public: a `git add .` on the cluster must not pick them up.

## Features

### Reuse the caller output of an earlier run
Run the pipeline again on the output of an earlier run without calling variants again, to try
other settings for merging, PoN, UNMET, MNV calling, annotation or post-filtering. Variant
calling is by far the most expensive step.

Points for the design:
- The raw calling depends on the BAM, reference, regions and `-Q` only. `-S`, `-D`, `-V` and `-M`
  are `bcftools filter` steps on the raw caller output (`lib/callers.sh`), so those can be
  changed too if the run starts from the raw caller files.
- Check that the earlier run used the same BAM, reference, regions and `-Q`, for example from
  a settings file written to `<output>/.circusnv/`, and refuse otherwise.
- Write to a new output directory and leave the earlier run untouched.
- `--cleanup` removes `isec/` and the `mnvcalling/` intermediates; check which files a reuse
  run needs, and keep those.

### MNVista-only run
Run only the MNV calling on an earlier run: build `mnvcalling/sites.vcf` from the stored
per-tool VCFs (or reuse it) and call MNVista, for example to try other MNVista settings or
`-B/--mnv-prior`. Shares most of its design with the reuse feature above.

## Iteration: biologically sound process

The restructuring kept the behaviour exactly as it was. In this iteration the process itself
is reviewed: each step and threshold should be a logical, biologically motivated choice.
Changes in results are allowed when they make the pipeline better, as separate commits with
the reasoning in the commit message and a Nibbler comparison that shows the effect.

Starting points from `docs/known-issues.md`: K1 strand-bias default, K2 UNMET threshold,
K5 matched-control filtering, K6 PoN percentage, D8 Mutect2 FILTER column, D9 shared PoN work
directory, D11 two annotation steps. Also the order of the steps (for example where PoN, UNMET
and post-filtering happen) and whether every caller-specific filter is still needed.

## Iteration: speed and resources

Make the pipeline as fast as possible, with sensible use of CPUs and memory. Speed matters
most. Every change is checked against the baseline: same results unless a change is meant to
alter them.

Starting points:
- D13: each caller gets the full thread count, so a sample uses about four times the requested
  CPUs. Divide the threads, or measure which caller is the bottleneck first.
- Run samples in parallel, or as a Slurm job array, instead of one after the other. Follow-up
  samples need the MNV results of their baseline (`-X`), so those have to wait for it.
- D11: the first OpenCRAVAT run uses all annotators while post-filtering needs only a few.
- Measure first: log the duration of every step per sample, so the slowest steps are known.
- Skip work that is not needed for the chosen options, and avoid writing and reading the same
  intermediate files more than once.
