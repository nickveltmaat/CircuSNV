# Known issues and open decisions

Each item needs a decision about the intended behaviour. Until it is decided, the pipeline
keeps the behaviour of the original code (the baseline), and the documentation describes that
behaviour. Items marked "decided" are done.

## Differences between the original README and the code

### K1. Default strand-bias threshold (`-S`)
- Code: 0.9 (the forward/reverse ratio must lie between 0.9 and 1.11).
- Original usage text and README: 0.1 (between 0.1 and 10).
- Now: the help text and README show the code default, 0.9.
- Hint: the `Strand Bias:` line in the settings block of old Slurm `.out` files shows which
  value validated runs used.
- Status: to be decided by Nick.

### K2. UNMET threshold
- Code: variants with an UNMET score of 0.90 or higher are removed.
- Original README: scores above 0.93 are removed.
- Now: README shows 0.90; the value is the setting `UNMET_THRESHOLD`.
- Status: to be decided by Nick.

### K3. Output directory (`-O`) and relative paths
- Decided: `-O` stays optional (default `./output`); relative paths resolve against the
  directory the pipeline is started from, now that the hard-coded `cd` is gone.

### K4. Python dependencies
- `requirements.txt` lists pandas, numpy, openpyxl, open-cravat (and scipy for
  `matched_filtering.py`), without versions.
- Open: pin the versions from `pip freeze` of the environment that produced the reference
  results.
- Status: waiting for that output.

### K5. Matched control (`-T`)
- Code: adds the reference/alternative read counts of the matched control to each variant;
  it does not filter. `scripts/matched_filtering.py` (Fisher test, VAF difference and ratio)
  is not called.
- Original README: described filtering with `matched_filtering.py`.
- Options: keep metrics only (README now says so); call `matched_filtering.py` again
  (behaviour change); or remove the script.
- Status: to be decided by Nick.

### K6. Percentage for single-tool PoN variants
- Code: 20% of the normals (setting `PON_MIN_NORMALS_PCT`). An old comment said 50%; the
  comment is gone.
- Status: to be decided by Nick (20 is used).

### K7. Directory paths with or without a trailing slash
- Decided and done: trailing slashes are removed once from `-I`, `-O` and `-P`, and paths are
  built as `dir/file` everywhere. SiNVICT still gets its directory arguments with a trailing
  slash, as before.

## Decisions found during the restructuring

### D8. Mutect2 FILTER column is ignored
- `FilterMutectCalls` runs, but the later `bcftools filter` steps do not look at the FILTER
  column, so variants that Mutect2 flags (for example as `weak_evidence`) are kept.
- Options: keep as is; or keep only `PASS` variants (behaviour change).
- Status: to be decided by Nick.

### D9. Panel-of-normals work directory is shared between runs
- The PoN is built in `PoN/` in the installation directory. The step that merges the normals
  lists every `PoN/*/` directory, so normals from an earlier run with other normals are
  included.
- Options: clear `PoN/` at the start of a PoN run, or build the PoN under the output directory
  (behaviour change when stale normals exist).
- Status: to be decided by Nick. Until then: remove `PoN/` before a run with a different set
  of normals.

### D10. Contig order in VCFs written by `isec_to_vcf.py`
- The `##contig` header lines come from a Python `set`, so their order can differ between runs.
  The records are not affected.
- Option: sort the contigs (changes only the header order).
- Status: to be decided by Nick.

### D11. Two annotation steps
- The first OpenCRAVAT run annotates with all annotators, although post-filtering only needs
  Sequence Ontology, 1000 Genomes/gnomAD frequencies and the repeat class.
- Option: annotate the first time with `thousandgenomes gnomad3 repeat` only. The final list
  stays the same, but `merged_vcfs_annotated.xlsx` and `removed_variants.xlsx` lose the
  columns of the other annotators.
- Status: to be decided by Nick (the two steps stay; see CLAUDE.md).

### D12. Matched control with a single BAM (`-T true -I file.bam`)
- Control BAMs are looked up in the input directory, so with a single BAM they are never
  found. The pre-run check now reports this as an error.
- Option: look up controls next to the single BAM, or add an option for a control directory.
- Status: to be decided by Nick.

### D13. Threads
- Each of the four callers gets the full thread count (`SLURM_CPUS_PER_TASK`), and SiNVICT
  starts as many bam-readcount processes, so a sample uses about four times the requested
  CPUs.
- Option: divide the threads over the callers. This may change the speed, and for some tools
  the order of records; it needs a test.
- Status: to be decided by Nick.

## To watch in the first test runs

- **Exit codes that were never checked.** With strict error handling, a non-zero exit of
  OpenCRAVAT (`oc run`), MNVista or a Python helper now stops the sample. If one of them
  returns non-zero although it worked, that call needs an exception.
- **SiNVICT input directory.** SiNVICT reads the bam-readcount directory, which holds both the
  per-chunk files and the merged file. Whether SiNVICT reads every file in it (and so sees
  each position twice) cannot be told from the code. The old and new version behave the same
  here. If it does, `<sample>_SiNVICT_RAW.vcf` would list calls twice; the normalised
  `<sample>_SiNVICT.vcf.gz` is deduplicated (`bcftools norm -d all`).
- **Flow diagram.** `flowchart_pipeline2.drawio.png` predates the restructuring and may not
  match every detail; it can only be updated in draw.io.
- **`raw-data/`.** Each sample gets an empty `raw-data/` directory, as before; the per-tool
  VCFs stay in the sample directory.
