# Conventions

Naming, message style and rules for pull requests in this repository.

## Contents

- Branch names
- Commit messages
- Pull requests
- Merge strategy
- Tools and quirks on Nick's machines
- Files that must never be committed

## Branch names

Format: `<type>/<short-kebab-description>`, lowercase, no spaces.

| Type | Use for |
| --- | --- |
| `feat` | new functionality |
| `fix` | a bug fix |
| `refactor` | restructuring with unchanged output |
| `docs` | documentation only |
| `config` | configuration, settings, skills |
| `test` | tests, test data, baseline tooling |
| `setup` | project or environment setup |
| `chore` | everything else that is housekeeping |

Existing examples: `setup/claude-code`, `docs/annotation-rationale`, `config/claude-permissions`.

## Commit messages

- English, imperative mood ("Add ...", "Fix ...", "Document ..."), subject at most 72 characters, no trailing period.
- Optional body after a blank line, explaining why the change is made rather than repeating the diff.
- One logical change per commit. A refactor and a behaviour change go in separate commits: that is what makes a baseline comparison meaningful, because a refactor must leave the output unchanged and a behaviour change must show up in it.
- The app adds a `Co-Authored-By` trailer by default. Leave it in place unless Nick asks otherwise.

## Pull requests

- Base branch `main`. One purpose per pull request.
- Title in the style of a commit subject.
- Description in English, based on `assets/pr-template.md`. Fill every section, remove the comments, and tick only what is true.
- State clearly whether the pipeline output can change ("Behaviour change" section). Nick relies on this to decide whether a Nibbler test is needed.
- Do not merge in the same step that opens the pull request.

## Merge strategy

Use a merge commit: `gh pr merge <N> --merge --delete-branch`. Do not squash or rebase unless Nick asks.

Reason: the history so far uses merge commits, and with a merge commit `git branch -d` works afterwards because the branch is part of `main`. A squash or rebase makes git consider the local branch unmerged.

## Tools and quirks on Nick's machines

- `git` runs through Git Bash. `gh` runs through PowerShell; it is not found in Git Bash.
- Multi-line text for `gh pr create`: write the body to a file outside the repository (for example in `$env:TEMP`) and pass it with `--body-file`. Alternatively use `--body` with a single-quoted PowerShell here-string. Do not create PR body files inside the repository.
- Line endings are LF, enforced by `.gitattributes` for scripts and Markdown. If a file has CRLF, fix it before committing; the preflight script reports this.
- Permission rules in `.claude/settings.json`: `git status`, `git diff` and `git log` run without asking; `git add`, `git commit`, `git push` and `gh pr create|merge|close` ask first; force-push, `reset --hard`, `clean`, `gh auth`, `gh api`, `gh secret` and `gh repo delete` are blocked. Other commands also ask under the default permission mode.

## Files that must never be committed

- Sequencing data: FASTQ, BAM, CRAM and their indexes.
- Pipeline output, logs and the `output/`, `results/` and `tmp/` directories.
- Anything with patient identifiers, sample sheets with identifiers, or pseudonymisation keys.
- Credentials, tokens, private keys.
- VCF files only if they are small synthetic test data that Nick agreed to commit.

If a file of this kind is already staged, unstage it (`git restore --staged <file>`) and tell Nick. If it was already pushed, stop and tell Nick right away: removing it from the history is a decision for Nick, not a quick fix.
