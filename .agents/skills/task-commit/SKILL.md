---
name: task-commit
description: Safely inspects working tree changes, checks for sensitive files, drafts a conventional commit message, presents a pre-commit report for user review, and executes the git commit upon approval. Handles both standard task commits and automated release workflows with Git tagging. Use whenever asked to commit changes, save progress to git, or create a release commit.
---

# Task Commit Skill

When triggered, act as a strict Git Commit and Release Manager. Follow this exact 5-step sequence:

---

## Step 1: Inspect (Silent)

1. **Detect Release Intent**:
   - Check if the user's prompt or task indicates a release workflow (e.g. *"commit for a new release"*, *"release commit"*, *"cut a release"*, *"commit v1.0.0"*, or `--release`).
   - If **Release Mode** is detected:
     a. **Run Version Check**:
        ```powershell
        D:\anaconda3\envs\pytransdat\python.exe tools/bump_version.py --check
        ```
        - If versions are mismatched or invalid, alert the user immediately and report the discrepancy. Offer to synchronize them using `D:\anaconda3\envs\pytransdat\python.exe tools/bump_version.py <version>` before proceeding.
     b. **Read Current Release Version**:
        - Retrieve the validated version `V` (e.g. `1.0.0-beta`, `1.0.0`).
     c. **Pre-flight Quality Gate Reminder**:
        - Verify that core geodetic precision (`regression-verify` or `pytest`) and frontend tests (`vitest`) have passed. Never commit a release if regressions or broken tests are outstanding.

2. **Check Working Tree Status**:
   - Run `git status -u` (or `git status --porcelain`). *(Note: On Windows sandbox environments where `.git` is isolated, run git inspection and commit commands with `BypassSandbox: true`)*.
   - If the working tree is completely clean (no modified, deleted, or untracked files), inform the user immediately:
     > *"Working tree is clean. There are no changes to commit."*
     and stop execution.

3. **Check Diffs & Existing Stage**:
   - Run `git diff` to inspect unstaged modifications.
   - Run `git diff --cached` (or `git diff --staged`) to check if any files were already staged prior to this invocation.

4. **Sensitive File & Cache Scan**:
   - Scan untracked and modified files for sensitive content or build artifacts:
     - Environment files (`.env`, `.env.*`)
     - Credentials, private keys, or tokens (`*.pem`, `*.key`, `id_rsa*`, secrets)
     - Caches and binary artifacts (`__pycache__`, `*.pyc`, `.pytest_cache/`, `*.log`, `Thumbs.db`, `.DS_Store`)
   - Flag any sensitive or cache files as **Excluded/Ignored**. Never propose staging them. If appropriate, recommend adding them to `.gitignore`.

---

## Step 2: Generate the Pre-Commit Report

Do **not** stage (`git add`) or commit anything yet. Present a clear, structured report to the user formatted as follows:

### Report Structure:
1. **Mode**: Standard Task Commit OR **Release Commit (`v<version>`)**
2. **Files to Stage**:
   - A bulleted list of the exact file paths you propose staging.
   - For each file, provide a concise 1-sentence summary of why it is being included based on the diff.
3. **Excluded / Ignored Files** *(omit if none)*:
   - List any untracked or cache files that are intentionally left unstaged, with a brief explanation.
4. **Proposed Commit Message**:
   - Formatted following the **Conventional Commits** standard:
     - **For Standard Commits**: `<type>(<optional-scope>): <concise description>`
       - Types: `feat`, `fix`, `refactor`, `test`, `docs`, `perf`, `chore`.
       - Optional scope: e.g., `test(trans_ro)`, `refactor(grid)`, `feat(helmert)`.
     - **For Release Commits**: `chore(release): v<version>`
       - Body: Summarize key features, bug fixes, breaking changes, or milestone tasks included in this release.
   - If working on a specific task file (e.g., in `.agents/tasks/`), reference the task identifier if relevant (e.g., `Task 06`).
5. **Git Tag to Create** *(Release Mode only)*:
   - Note the proposed annotated tag: `v<version>` (e.g. `v1.0.0-beta`, `v1.0.0`).

---

## Step 3: Pause for Review (Interactive Buttons)

After presenting the pre-commit report in the chat, **do NOT proceed to Step 4**.
Call the `ask_question` tool so the Antigravity IDE renders interactive clickable buttons for the user to approve, edit, or cancel:

### A. For Standard Commits:
- **Question**: `"Review the pre-commit report above. How would you like to proceed?"`
- **Options**:
  1. `"(Recommended) Approve and execute the commit as proposed"`
  2. `"Edit the commit message before committing"`
  3. `"Cancel the commit"`
- **is_multi_select**: `false`

### B. For Release Commits:
- **Question**: `"Review the release pre-commit report above. How would you like to proceed?"`
- **Options**:
  1. `"(Recommended) Approve release commit and create Git tag v<version>"`
  2. `"Approve release commit only (do not create Git tag)"`
  3. `"Edit the commit message before committing"`
  4. `"Cancel the release commit"`
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

3. **Release Git Tagging (Release Mode only, upon user approval)**:
   - If the user approved creating the Git tag:
     ```powershell
     git tag -a v<version> -m "PyTransDatRO v<version>"
     ```
   - Verify that the tag was created:
     ```powershell
     git tag -l -n1 v<version>
     ```

4. **Pushing Policy**:
   - **NEVER** run `git push`. Leave remote pushing strictly to the user.

---

## Step 5: Post-Commit Verification

After committing:
1. Run `git log -1 --stat` to retrieve the created commit details.
2. If in Release Mode and tagged, run `git tag -l -n1 v<version>` to confirm the tag.
3. Run `git status` to verify the state of the remaining working tree.
4. Present a clear confirmation to the user with:
   - The commit hash and subject.
   - The created Git tag (if applicable).
   - Any remaining unstaged files.
   - **For Release Mode**: Output the exact command for the user to push the release and tags to the remote repository:
     ```powershell
     git push origin main --tags
     ```
