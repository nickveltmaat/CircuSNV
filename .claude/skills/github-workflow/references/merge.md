# Merge a pull request

Merging is the last step and it is irreversible in practice: the change lands on `main`. Do it only on Nick's go-ahead, one command per tool call.

## Preconditions

- Nick has said that the Nibbler test passed, or has confirmed that the change cannot affect pipeline behaviour (documentation or configuration only). If this is not clear, ask. Do not infer it.
- The pull request number is known. If Nick did not give one, run `gh pr list` and ask which one.

## Steps

Use PowerShell for `gh` commands and Git Bash for `git` commands.

1. `gh pr view <N>`: check that the state is open, the base is `main`, and the title and branch are the ones Nick means.
2. Tell Nick in Dutch what is about to land on `main` (title, number, files changed) and wait for a clear yes.
3. `gh pr merge <N> --merge --delete-branch`. Use a merge commit, not squash or rebase (see `references/conventions.md`).
4. Verify locally, one command each:
   - `git status`
   - `git log --oneline --graph -5`
   - `git branch -a`
5. Clean up what is still left:
   - If the working copy is not on `main` or `main` is behind: `git switch main`, then `git pull --prune`.
   - If the local branch still exists: `git branch -d <branch>`. Never `-D`. If `-d` refuses because the branch is not merged, stop and report: something is not as expected.
6. Report to Nick in Dutch: the PR is merged, `main` is at `<short hash>`, the branch is gone locally and remotely, and what the next step is.

## When something goes wrong

- Merge conflicts or a PR that cannot be merged: stop and report the message. Do not resolve conflicts without Nick.
- `gh` is not found: you are in the wrong shell. Run it in PowerShell.
- Authentication problems: stop. Logging in is something Nick does himself with `gh auth login`.
