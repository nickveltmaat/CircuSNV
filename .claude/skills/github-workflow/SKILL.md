---
name: github-workflow
description: Git and GitHub workflow for the CircuSNV repository (nickveltmaat/CircuSNV) - create a branch, run preflight checks, commit, push, open a pull request, hand a branch over for testing on the Nibbler HPC cluster, and merge. Use this whenever the user asks to commit, push, make or merge a pull request, create or clean up a branch, "ship" a change, or prepare a change for a Nibbler test, even if they do not mention GitHub explicitly.
argument-hint: "[ship | merge <PR number> | status]"
---

# GitHub workflow for CircuSNV

This skill describes how a change travels from a local branch to `main` in this repository. Nick approves every step in the app, and the test on the Nibbler cluster is done by Nick himself.

## Repository facts

- Repository: `nickveltmaat/CircuSNV` (https://github.com/nickveltmaat/CircuSNV), public, default branch `main`, HTTPS remote.
- The repository is public, so everything committed is readable by anyone. Patient data, sequencing data, credentials and tokens must never end up in a commit. `.gitignore` is a safety net, not a plan: stage files by name.
- On Nick's Windows machines `git` runs through Git Bash and `gh` runs through PowerShell (`gh` is not on the Git Bash PATH). Use the right tool for each.

## How to pick a procedure

| The user wants to... | Read |
| --- | --- |
| commit, push, open a PR, "ship" a change (`/github-workflow ship`) | `references/ship.md` |
| test a branch on Nibbler | `references/cluster-handoff.md` |
| merge a PR (`/github-workflow merge 12`) | `references/merge.md` |
| naming, commit messages, PR rules, merge strategy | `references/conventions.md` |
| `status` | run `git status` and `git log --oneline --graph -5`, then `gh pr status`, and summarise |

Read only the file that matches the request. Each file is self-contained.

## Rules that always apply, and why

- **Never work on `main`.** The pull request is the review moment, and `main` has to stay a known-good baseline for the pipeline. If the current branch is `main`, stop and propose a branch name.
- **One command per tool call.** Do not chain git or gh commands with `&&`, `;` or pipes. Nick approves each git/gh action in the app, and a chained command hides what is being approved.
- **Stage files by name.** Never `git add -A`, `git add .` or wildcards: they pull in data and output by accident.
- **Do not change permissions.** Never edit `.claude/settings.json` or `.claude/settings.local.json` on your own, never suggest "always allow", and never run `gh auth`, `gh api`, `gh secret` or `gh repo delete`. Those are Nick's to do.
- **Never merge before the test result is in.** Merge only after Nick says the Nibbler test passed, or after Nick confirms the change cannot affect pipeline behaviour (documentation or configuration only).
- **If something is refused or fails, stop.** Report the command and the error. Do not look for a workaround: a refusal means a rule did its job, or something is wrong that Nick should see.
- **No history rewriting.** No force-push, `reset --hard`, `clean` or rebase, and no amending of commits that are already pushed, unless Nick explicitly asks.

## Language

Explain things to Nick in Dutch. Branch names, commit messages, PR titles and PR descriptions are in English because the repository is public.
