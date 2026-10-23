# Ship: from local change to pull request

Take a local change to an open pull request. Do the steps in order, one command per tool call, and stop at the end of step 9.

## Contents

- Step 1: Orient
- Step 2: Preflight
- Step 3: Review the diff
- Step 4 to 7: Stage, commit, push
- Step 8: Open the pull request
- Step 9: Report and stop

## Step 1: Orient

Run `git status`, then `git branch --show-current`.

- On `main`: stop. Propose a branch name following `references/conventions.md` and wait for Nick to agree. Then run `git switch -c <name>`. Uncommitted changes move along to the new branch.
- On another branch: continue.
- Changes that clearly belong to different purposes: ask Nick whether to split them into separate branches or commits before going on.

## Step 2: Preflight

Run `bash .claude/skills/github-workflow/scripts/preflight.sh`. It only reads; it changes nothing.

It checks the branch, data files, secrets, line endings, bash and Python syntax, JSON validity and file size.

- `FAIL` lines block the commit. Fix the cause (or ask Nick if the fix is unclear) and run the script again.
- `WARN` lines do not block. Mention them to Nick.
- `NOTE` lines are reminders, for instance that pipeline code changed.

VCF files fail by default, because the repository is public. Only if Nick explicitly agrees that a VCF is small synthetic test data, run the script as `PREFLIGHT_ALLOW_VCF=1 bash .claude/skills/github-workflow/scripts/preflight.sh`.

## Step 3: Review the diff

Run `git diff --stat`, then `git diff`. For files that are new and not yet tracked, `git status --short -uall` lists them and the files can be read directly.

Tell Nick, in Dutch, in a few sentences:
- what changed and why;
- whether pipeline output or defaults can change (a behaviour change), yes or no;
- if the diff mixes a refactor with a behaviour change, propose splitting it into two commits.

## Step 4 to 7: Stage, commit, push

4. `git add <file> <file>`: only files that belong to this change, by name.
5. `git diff --cached --stat`: check that exactly these files are staged.
6. `git commit -m "<subject>"`, with an extra `-m "<body>"` when the reason needs explaining. Follow `references/conventions.md`.
7. `git push -u origin <branch>`.

If the push is refused or asks for authentication, stop and report. Do not retry in a different way.

## Step 8: Open the pull request

Use PowerShell, because `gh` is not available in Git Bash.

1. Fill in `assets/pr-template.md` for this change and save the result outside the repository, for example as `$env:TEMP\pr-body.md`. Remove the comment lines.
2. Run `gh pr create --base main --head <branch> --title "<title>" --body-file "$env:TEMP\pr-body.md"`.

If writing a temporary file is not possible, pass the text with `--body` and a single-quoted PowerShell here-string.

## Step 9: Report and stop

Tell Nick, in Dutch:
- the pull request link;
- the answer to the question "can pipeline behaviour change?";
- if it can: that a Nibbler test is needed before merging, followed by the commands from `references/cluster-handoff.md` for this branch;
- if it cannot (documentation or configuration only): that no cluster test is needed and the pull request can be reviewed and merged with `/github-workflow merge <number>` once Nick agrees.

Do not merge. Merging is a separate step in `references/merge.md`.
