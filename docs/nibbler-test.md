# Testing a CircuSNV branch on Nibbler

This procedure runs a fixed baseline version and a candidate version of CircuSNV on the
same small, anonymised input and compares their output record by record. Use it for every
phase of the `refactor/circusnv-v2` work: a refactor must give identical output, and an
approved behaviour change must show only the expected difference.

No test data is stored in this repository. Inputs stay on the cluster.

## Baseline

The baseline is commit `e88243c8361e53b107400ce2e662e6e767189b1e` on `main`. Its pipeline
code is unchanged since commit `9dc46f6`; later commits only added documentation and
configuration.

## Placeholders

Replace these everywhere below. They are not stored in the repository on purpose.

| Placeholder | Meaning |
| --- | --- |
| `<group>`, `<tmp>` | your group and its `tmp*` area: Slurm only accepts jobs whose working directory, `--output` and `--error` are under `/groups/<group>/<tmp>/` |
| `<test_dir>` | a directory under that area for these tests |
| `<install_dir>` | your existing, working CircuSNV installation (with `env`, `MNVista` and resources) |
| `<python_module>` | the module that provides the Python used by `<install_dir>/env`, if one is needed |
| test input | a small, anonymised BAM (or directory of BAMs), reference FASTA and BED, and optionally normals and a crosslink table |

## 1. One-time setup

On a Nibbler user interface server, one line at a time:

```
module load git
mkdir -p /groups/<group>/<tmp>/<test_dir>/logs /groups/<group>/<tmp>/<test_dir>/runs
cd /groups/<group>/<tmp>/<test_dir>
git clone https://github.com/nickveltmaat/CircuSNV.git
cd CircuSNV
git worktree add --detach ../circusnv-baseline e88243c8361e53b107400ce2e662e6e767189b1e
git worktree add --detach ../circusnv-candidate origin/refactor/circusnv-v2
```

Copy the sbatch template and fill in its `<placeholders>`:

```
cp tests/nibbler/run_circusnv.sbatch.template ../run_circusnv.sbatch
nano ../run_circusnv.sbatch
```

## 2. Prepare the worktrees

Both worktrees need the resources that are not in the repository. Link them from your
installation; skip a line if you do not use that feature:

```
cd /groups/<group>/<tmp>/<test_dir>
for wt in circusnv-baseline circusnv-candidate; do ln -s <install_dir>/env $wt/env; done
for wt in circusnv-baseline circusnv-candidate; do ln -s <install_dir>/MNVista $wt/MNVista; done
for wt in circusnv-baseline circusnv-candidate; do ln -s <install_dir>/MNV_PoN_No_overlap_new.txt $wt/MNV_PoN_No_overlap_new.txt; done
for wt in circusnv-baseline circusnv-candidate; do ln -s <install_dir>/hg38_genome.UNMET_bin10bp.bed $wt/hg38_genome.UNMET_bin10bp.bed; done
```

The other tools (LoFreq, SiNVICT, bam-readcount, VarDict) come from your `PATH`, which the
job inherits through `--get-user-env`.

### Patch the hard-coded `cd` (baseline always; candidate until phase 2)

The old `CircuSNV.sh` changes to a fixed directory, so it would run the scripts of that
directory instead of the worktree's. In each worktree that still has this line, replace it
locally (never commit this) with a change to the script's own directory:

```
sed -i 's|^cd /groups/.*$|cd "$(dirname "$(readlink -f "$0")")"|' circusnv-baseline/CircuSNV.sh
git -C circusnv-baseline diff --stat
```

The diff must show `1 file changed, 1 insertion(+), 1 deletion(-)`. Do the same for
`circusnv-candidate` as long as `grep -n '^cd /groups/' circusnv-candidate/CircuSNV.sh`
prints a line.

## 3. Choose the commit to test

Every status update and the pull request list the commit hash of each finished phase.

```
cd /groups/<group>/<tmp>/<test_dir>/CircuSNV
git fetch origin
git -C ../circusnv-candidate checkout --detach <phase_commit>
git -C ../circusnv-candidate log --oneline -3
```

If the candidate worktree has the local `cd` patch, run `git -C ../circusnv-candidate checkout -- CircuSNV.sh`
before the checkout and apply the patch again afterwards if it is still needed.

## 4. Run both versions

Use absolute paths for every input, and exactly the same options for both runs. Do not pass
`-O`: the output directory is the second argument. Each run needs a new, empty output directory.

```
cd /groups/<group>/<tmp>/<test_dir>
rm -rf circusnv-baseline/PoN circusnv-candidate/PoN
sbatch run_circusnv.sbatch /groups/<group>/<tmp>/<test_dir>/circusnv-baseline /groups/<group>/<tmp>/<test_dir>/runs/baseline-1 -I <bam_or_dir> -R <reference.fasta> -L <regions.bed>
sbatch run_circusnv.sbatch /groups/<group>/<tmp>/<test_dir>/circusnv-candidate /groups/<group>/<tmp>/<test_dir>/runs/candidate-1 -I <bam_or_dir> -R <reference.fasta> -L <regions.bed>
```

Removing `PoN` matters when you test with `-P`: the old pipeline merges every normal it
finds under `./PoN`, including those from earlier runs.

Suggested runs, from small to complete: first the minimal command above; then with `-P`;
then with the options you use in practice (`-H`, `-X`, `-T`, `-U`).

Follow the jobs:

```
squeue -u $USER
sacct -j <job_id> --format=JobID,State,ExitCode,Elapsed,MaxRSS,AllocCPUS
```

## 5. Compare

The compare script needs `openpyxl` for the Excel files; the pipeline venv has it. For
small test output this is light enough for the user interface server; for large output
use `srun`.

```
module load <python_module>
cd /groups/<group>/<tmp>/<test_dir>
circusnv-candidate/env/bin/python3 circusnv-candidate/tests/nibbler/compare_outputs.py runs/baseline-1 runs/candidate-1
echo "exit status: $?"
```

Exit status 0 means identical, 1 means differences, 2 means a file could not be compared.

If the candidate is still at a commit without `tests/nibbler/`, run the script from
`CircuSNV/tests/nibbler/` instead (after `git -C CircuSNV switch refactor/circusnv-v2`).

### Measure the noise first (phase 1)

In phase 1 the candidate has the same pipeline code as the baseline. Compare those two runs
first: any difference is non-determinism of the pipeline or its tools, not of the refactor.
Report it, so later phases can be judged against it.

## 6. What to paste back

Safe to paste:
- the full report of `compare_outputs.py` as printed by default (sample directories are
  masked as S1, S2, ...);
- the `sacct` line of each job;
- the `Exit status` line and the last 20 to 30 lines of the job's `.out` and `.err` files,
  after checking them for sample names, patient numbers and project paths.

Never paste: VCF or Excel content, sample names or tables with identifiers, the output of
`--no-mask`, `--show-paths` or `--show-records`, or the mapping file of `--mapping-file`.

## 7. Clean up

```
cd /groups/<group>/<tmp>/<test_dir>/CircuSNV
git worktree remove --force ../circusnv-candidate
git worktree remove --force ../circusnv-baseline
```

`--force` is needed because of the local `cd` patch and the resource links.
