---
name: task-commit
description: Safely inspects working tree changes, checks for sensitive files, drafts a conventional commit message, presents a pre-commit report for user review, and executes the git commit upon approval. Use whenever asked to commit changes, save progress to git, or create a task commit.
---

# Task Commit Skill

When triggered, act as a strict Git Commit Manager. Follow this exact 5-step sequence:

---

## Step 1: Inspect (Silent)

1. **Check Working Tree Status**:
   - Run `git status -u` (or `git status --porcelain`).
   - If the working tree is completely clean (no modified, deleted, or untracked files), inform the user immediately:
     > *"Working tree is clean. There are no changes to commit."*
     and stop execution.

2. **Check Diffs & Existing Stage**:
   - Run `git diff` to inspect unstaged modifications.
   - Run `git diff --cached` (or `git diff --staged`) to check if any files were already staged prior to this invocation.

3. **Sensitive File & Cache Scan**:
   - Scan untracked and modified files for sensitive content or build artifacts:
     - Environment files (`.env`, `.env.*`)
     - Credentials, private keys, or tokens (`*.pem`, `*.key`, `id_rsa*`, secrets)
     - Caches and binary artifacts (`__pycache__`, `*.pyc`, `.pytest_cache/`, `*.log`, `Thumbs.db`, `.DS_Store`)
   - Flag any sensitive or cache files as **Excluded/Ignored**. Never propose staging them. If appropriate, recommend adding them to `.gitignore`.

---

## Step 2: Generate the Pre-Commit Report

Do **not** stage (`git add`) or commit anything yet. Present a clear, structured report to the user formatted as follows:

### Report Structure:
1. **Files to Stage**:
   - A bulleted list of the exact file paths you propose staging.
   - For each file, provide a concise 1-sentence summary of why it is being included based on the diff.
2. **Excluded / Ignored Files** *(omit if none)*:
   - List any untracked or cache files that are intentionally left unstaged, with a brief explanation.
3. **Proposed Commit Message**:
   - Formatted following the **Conventional Commits** standard:
     `<type>(<optional-scope>): <concise description>`
     - Types: `feat`, `fix`, `refactor`, `test`, `docs`, `perf`, `chore`.
     - Optional scope: e.g., `test(trans_ro)`, `refactor(grid)`, `feat(helmert)`.
   - Include a body paragraph if the changes are complex or span multiple functional areas.
   - If working on a specific task file (e.g., in `.agents/tasks/`), reference the task identifier if relevant (e.g., `Task 06`).

---

## Step 3: Pause for Review (Interactive Buttons)

After presenting the pre-commit report in the chat, **do NOT proceed to Step 4**.
Call the `ask_question` tool so the Antigravity IDE renders interactive clickable buttons for the user to approve, edit, or cancel:

- **Question**: `"Review the pre-commit report above. How would you like to proceed?"`
- **Options**:
  1. `"(Recommended) Approve and execute the commit as proposed"`
  2. `"Edit the commit message before committing"`
  3. `"Cancel the commit"`
- **is_multi_select**: `false`

Handle the response:
- If the user selects **Approve**: Proceed to Step 4.
- If the user selects **Edit**: Ask the user for their desired commit message or suggest revisions, then prompt again.
- If the user selects **Cancel**: Terminate the workflow without making any commits.

Wait for explicit user approval via the `ask_question` tool before proceeding to Step 4. **Never proceed to Step 4 in the same turn.**


---

## Step 4: Execute (Upon Explicit Approval)

Once the user approves:

1. **Targeted Staging**:
   - Run `git add <file1> <file2> ...` targeting **only** the specific approved files.
   - **CRITICAL**: Never run blanket commands like `git add .` or `git add -A`.

2. **Safe Commit Execution (PowerShell-Safe)**:
   - On Windows/PowerShell, avoid single multi-line strings in `-m "..."` that can break on quotes or newlines.
   - Use either:
     - Multiple `-m` arguments:
       ```powershell
       git commit -m "<subject>" -m "<body>"
       ```
     - Or write the approved commit message to a temporary text file in the scratch directory and commit with:
       ```powershell
       git commit -F <path_to_temp_commit_msg>
       ```

3. **Pushing Policy**:
   - **NEVER** run `git push`. Leave remote pushing strictly to the user.

---

## Step 5: Post-Commit Verification

After committing:
1. Run `git log -1 --stat` to retrieve the created commit details.
2. Run `git status` to verify the state of the remaining working tree.
3. Present a brief confirmation to the user with the commit hash, commit message subject, and whether any unstaged files remain.
