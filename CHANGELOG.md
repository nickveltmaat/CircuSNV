# Changelog

## 2.0.0-dev (unreleased)

A restructured version with the same scientific behaviour: the callers, filters, thresholds
and output files are unchanged unless listed under "Behaviour changes" or "Fixes".

### Behaviour changes
- The MNV panel of normals is off by default and is used only with `--mnv-pon FILE`. Before,
  MNVista always got `MNV_PoN_No_overlap_new.txt` from the installation directory whenever
  `-P` was given; now `-P` is only the PoN for the SNV calls. To get the old MNV results with
  `-P`, add `--mnv-pon <install>/MNV_PoN_No_overlap_new.txt`.
- The per-tool calls (main path and MNV input) are restricted to the regions of `-L`, after
  normalisation: a call is kept when its REF span overlaps a region. LoFreq called on whole
  reads that touch a region, so it reported positions up to a read length outside the panel,
  and VarDict reported some calls within a few bases of a region. These calls are gone from
  the per-tool VCFs, the merged lists, the MNV input and the final list. The run log shows
  per sample and tool how many calls were removed. The normals of the panel of normals are
  restricted the same way.
- With a panel of normals, an indel is also removed when the PoN has an indel with other
  alleles at the same position (`--pon-indel-match position`, the default). The PoN filter
  compared exact alleles only, so a 1-bp insertion in a short G/C stretch passed when the
  normals showed the opposite 1-bp deletion at the same anchor base (polymerase slippage).
  SNVs and the MNV input are not affected. The removed indels are listed per sample in
  `intermediate-files/merged_vcfs_PoN_indel_position_removed.txt`;
  `--pon-indel-match allele` gives the old behaviour.

### Structure
- `CircuSNV.sh` is a short entry point; the code lives in `lib/` (one module per stage) and
  the Python helpers in `scripts/`.
- No hard-coded `cd`: the installation directory is found automatically (also under
  `sbatch`). Relative paths, including the default `./output`, are taken from the directory
  the pipeline is started from.
- Module versions, resource paths and fixed values are settings with the old values as
  defaults, overridable with `--config` (see `config/circusnv.example.conf`).
- Variables are quoted throughout; MNVista is called without `eval`.

### Usability
- Long options next to every short option, `--help` with the real defaults, `--version`.
- All options are validated before starting, and every problem is reported at once.
- Pre-run checks for inputs, the Python environment, MNVista, modules and commands
  (`--skip-checks` to turn them off).
- `--dry-run` shows the settings, samples and steps without running anything.
- Timestamped output and a log file per run in `<output>/.circusnv/`.
- `--resume` skips samples that finished in an earlier run; `--cleanup` removes
  intermediate files after each sample; `--threads` sets the thread count.
- `--mnv-pon` sets the MNV panel of normals (see "Behaviour changes").
- `--pon-indel-match position|allele` chooses how the PoN filter matches indels (see
  "Behaviour changes"); the settings block shows the value.
- `--mnv-flank-regions PATTERN` and `--mnv-flank-bp INT` (default 200): with `-H true`, the
  regions whose last BED column matches the pattern are extended for calling and the MNV
  input, so MNVista also sees variants in their flanks; the main results stay inside the
  panel. `--dry-run` and the run log show how many regions match.

### Robustness
- The pipeline runs with `set -euo pipefail`. The exit status of every caller is checked, for
  tumor samples and for the normals of the panel of normals (which now also get retries).
- A failing sample no longer stops the other samples; the run ends with a list of failed
  samples and a non-zero exit status.
- The Python helpers exit with status 1 on wrong arguments.

### Fixes
- The per-tool PoN lists could be read before they were written (race condition).
- VarDict's cleanup could delete intermediate files of LoFreq and SiNVICT while they ran.
- A BAM index could be taken for present when another sample's name started the same way.
- Retries of a caller failed on files left by the previous attempt.
- LoFreq's temporary BAM stayed behind when LoFreq failed.
- A missing MNV panel of normals (a `--mnv-pon` file that does not exist, with `-H true`)
  stops the run before it starts; before, MNVista got the missing file without a warning.
- `isec_to_vcf.py` is now the version that runs on the cluster. The version in git was older:
  it wrote the MNV input without AF, DP and MRD as INFO fields, which gave different MNV
  results. The contigs in its VCF headers are now sorted.

### Efficiency
- Normalised caller VCFs are written compressed directly.
- The four MNV inputs are normalised in parallel.

### Documentation and testing
- README rewritten; `requirements.txt`; `docs/known-issues.md` lists open decisions.
- `tests/nibbler/compare_outputs.py` compares the output of two runs record by record;
  `docs/nibbler-test.md` describes the test on the cluster.
