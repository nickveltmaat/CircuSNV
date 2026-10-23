# Handoff to the Nibbler cluster

Claude has no access to the Nibbler cluster. Nick runs the test there himself. Claude's job is to give him exact commands for the branch at hand and to interpret what he brings back. That keeps credentials and data on the cluster and out of this session.

## Contents

- When a cluster test is needed
- Commands to give Nick
- The test itself
- What Nick may paste back
- Rules for anything written for Nibbler
- After the result

## When a cluster test is needed

A test is needed when the change can alter pipeline behaviour: edits to `CircuSNV.sh` or any Python script it calls, changed defaults or parameters, changed tool versions. It is not needed for documentation, `CLAUDE.md`, skills or permission settings. Say explicitly which of the two applies.

## Commands to give Nick

Give these as a code block, with the real branch name filled in. Nick runs them on a Nibbler user interface server, one line at a time.

```
module load git
cd <path to Nick's clone of CircuSNV on Nibbler>
git fetch origin
git switch <branch>
git log --oneline -3
```

The last command should show the latest commit of the pull request at the top. If Nick has no clone on the cluster yet, the first time he needs `git clone https://github.com/nickveltmaat/CircuSNV.git` in a suitable directory. If the branch already exists on the cluster and new commits were pushed, he adds `git pull` after the switch.

`module load git` matters: the default git on the cluster can be too old for https access to GitHub.

## The test itself

The test procedure (which small dataset, which command, which baseline output to compare with) is not defined yet. If `docs/nibbler-test.md` exists in the repository, follow it and give Nick its commands. If it does not exist, say so plainly and do not invent a test command or a comparison. Defining the baseline is the next piece of work for this project.

## What Nick may paste back

Useful and safe to share:
- the exit status and the last 20 to 30 lines of the Slurm `.out` and `.err` files;
- counts and summaries (number of variants per stage, a diff summary against the baseline);
- the result of a comparison script.

Not to be pasted: raw sequencing data, complete VCFs of patient samples, file names or tables containing patient identifiers or patient numbers. The UMCG code of conduct for the cluster only allows pseudonymised or anonymised data and no identifiers. If Nick pastes something that looks like an identifier, say so and ask him to leave it out.

## Rules for anything written for Nibbler

These come from the cluster documentation and from errors seen earlier:
- Slurm jobs are rejected when the working directory, `--output` or `--error` is not under the group's `/groups/<group>/tmp*/` area.
- No resource-intensive computation on the user interface (login) server: use `sbatch`, or `srun` for interactive work.
- In job scripts use `#SBATCH --nodes=1`, `#SBATCH --export=NONE` and `#SBATCH --get-user-env`, load modules inside the script and finish with `module list` so the job log records the versions.
- Request resources that match what the job uses: the number of `--cpus-per-task` equal to the threads started.

## After the result

- Test passed: tell Nick the pull request can be merged and that `/github-workflow merge <number>` does it.
- Test failed: analyse the pasted output, fix on the same branch with a new commit, push, and give Nick the commands again (with `git pull` on the cluster).
- Output differs from the baseline: do not call it a failure or a success on your own. Describe the differences and let Nick decide whether the new behaviour is intended, then record that in the pull request.
