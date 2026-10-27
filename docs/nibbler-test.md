# Testing a CircuSNV branch on Nibbler

You need: one small, anonymised input that you ran before with the old pipeline, that run's
output directory (the reference) and its options (see the settings block in its
`CircuSNV.out`). Replace `<...>` with your own paths; `<work>` is a directory under your
group's `/groups/<group>/tmp*/` area, `<install>` is your current CircuSNV installation.

## Quick test

1. `module load git`
2. `cd <work>`
3. First time: `git clone https://github.com/nickveltmaat/CircuSNV.git circusnv-v2`
4. `cd circusnv-v2`
5. `git fetch origin` and `git switch refactor/circusnv-v2` (later updates: `git pull`)
6. Link the files that are not in git: `ln -s <install>/env env`, and when used also
   `ln -s <install>/MNVista MNVista`, `ln -s <install>/MNV_PoN_No_overlap_new.txt .`,
   `ln -s <install>/hg38_genome.UNMET_bin10bp.bed .` (or pass the MNV panel of normals with
   `--mnv-pon <install>/MNV_PoN_No_overlap_new.txt`)
7. Check without running: `bash CircuSNV.sh --dry-run <options of the reference run> -O <work>/new_output`
8. Run: `sbatch CircuSNV.sh <options of the reference run> -O <work>/new_output`
9. Follow it with `squeue -u $USER`; when it is gone, read `tail -n 30 CircuSNV.out CircuSNV.err`
10. Compare: `env/bin/python3 tests/nibbler/compare_outputs.py <reference_output> <work>/new_output`
11. Report back: the `Result:` line and the table of the comparison, the exit status, and the
    last lines of `CircuSNV.err`. Never paste sample names, patient numbers, VCF or Excel
    content, or the output of `--no-mask`/`--show-records`.

`Result: IDENTICAL` means the new version reproduces the reference. If `python3` is not found
in step 10, load the Python module that the env was built with first.

## Only if there are differences

Differences can come from the new code, or from things that changed since the reference run
(annotation databases, tool versions). To tell them apart, run the old version now on the same
input and compare the two fresh runs:

1. `cd <work>/circusnv-v2` and `git worktree add --detach ../circusnv-old e88243c8361e53b107400ce2e662e6e767189b1e`
2. In `../circusnv-old`, link `env`, `MNVista` and the other files as in step 6. That commit
   has an older `isec_to_vcf.py` than the cluster installation (different MNV results), so
   also copy it: `cp <install>/isec_to_vcf.py ../circusnv-old/`
3. The old script changes to a fixed directory; make it use the directory it is submitted
   from instead (local change, never committed):
   `sed -i 's|^cd /groups/.*$|cd "$SLURM_SUBMIT_DIR"|' ../circusnv-old/CircuSNV.sh`
   and check that `git -C ../circusnv-old diff --stat` shows 1 line changed.
4. Submit it from its own directory, with absolute paths for all inputs: `cd ../circusnv-old`
   and `sbatch CircuSNV.sh <options> -O <work>/old_output`
5. Compare: `../circusnv-v2/env/bin/python3 ../circusnv-v2/tests/nibbler/compare_outputs.py <work>/old_output <work>/new_output`
6. Clean up afterwards: `git -C <work>/circusnv-v2 worktree remove --force ../circusnv-old`

Report the same lines as above, for both comparisons.
