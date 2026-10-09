---
name: task-record
description: Creates or updates a task documentation file in .agents/tasks/ capturing the cumulative outcome, technical architecture, and verification of the current chat session. Use whenever asked to record, document, or update the current task or session results.
---

# Task Record Skill

## Purpose
This skill captures the cumulative work, architectural decisions, and verification results of the current chat session into a dedicated Markdown task file under `.agents/tasks/`. 

The produced task file serves as the definitive reference for any future agent or developer to understand:
- What was implemented or resolved.
- How architectural decisions were made and why.
- Which components/files were created or modified.
- How the implementation was verified and tested.

---

## Operating Principles

1. **Outcome-Oriented (The Final State)**:
   - Document the *final result* of the work achieved in the session.
   - Do **not** log conversational back-and-forth, intermediate bugs that were subsequently fixed, or iterative debate.
   - The file should read as an authoritative technical record written once upon completion.

2. **Session Grounding**:
   - Derive facts strictly from what occurred in the current chat:
     - Files modified or created (from edits/writes performed in the session).
     - Test runs and results (e.g., exact test names and pass rates like `70/70 passed`).
     - Architectural invariants discussed and maintained (e.g., standard library purity, isolation rules).

3. **Single File Per Session (Create vs. Update)**:
   - **First invocation in a session**:
     - Inspect `.agents/tasks/` to determine the highest existing two-digit numeric prefix (e.g., `09` -> next is `10`).
     - Pick a clear, descriptive lowercase snake_case slug summarizing the task (e.g., `10_feature_or_bugfix_description.md`).
     - Create the new task file.
   - **Subsequent invocations in the same session**:
     - Do **not** increment the prefix to create a new file.
     - Identify the task file already created for the current session's work and update it in-place to reflect the expanded or updated final outcome.

4. **No Meta Noise**:
   - Never reference the `task-record` skill itself, prompt instructions, or internal agent workflows in the task file.

---

## File Anatomy & Conventions

Task files in `.agents/tasks/` follow a clean, consistent structure while remaining flexible to the specific nature of the work:

```markdown
# <Clear, Descriptive Task Title>

## Objectives
- [x] <High-level objective 1 achieved in this session>
- [x] <High-level objective 2 achieved in this session>
- [ ] <Optional: explicitly deferred or planned follow-up if discussed>

## Implementation Details

1. **<Component / Focus Area 1> (`<path/to/file>`)**:
   - Detailed explanation of changes, classes, or functions introduced or modified.
   - Design rationale: why this approach was chosen.
   - Invariants respected (e.g., pure Python standard library, error handling, performance).

2. **<Component / Focus Area 2> (`<path/to/file>`)**:
   - Details of data structures, algorithms, or protocol handling.

## Verification & Testing
- **Test Modules**: `<tests/test_foo.py>`, `<tests/test_bar.py>`
- **Scenarios Covered**: Brief list of edge cases, benchmarks, or workflows verified.
- **Results**: Final test pass rate (e.g., `48/48 core passed, 22/22 API passed (70/70 total)`).
```

*Note*: If the task is of a distinct nature (e.g., pure benchmarking, research investigation, or refactoring cleanup), adapt the headings as needed to best communicate the outcome, but maintain clear technical depth.

---

## Workflow Steps

When triggered, execute the following steps:

### Step 1: Scan Existing Tasks & Determine Target Path
1. List `.agents/tasks/` to inspect all existing filenames (`01_...`, `02_...`, ..., `NN_...`).
2. Determine whether a task file was already created for the current session:
   - If yes: Select that file for an in-place update.
   - If no: Take the maximum prefix number, increment by 1 (formatted as two digits with leading zero, e.g., `10`), and format the file path as `.agents/tasks/<NN>_<descriptive_slug>.md`.

### Step 2: Synthesize Session Work
Gather the following elements from the chat conversation:
- **Core Goal**: What was the overarching objective?
- **Files Modified / Created**: Identify all workspace files touched during the session.
- **Architectural & Design Choices**: Extract key design patterns, invariants preserved, and rationale agreed upon.
- **Verification**: Note which test commands were run and their outcome/pass rates.

### Step 3: Write or Update the Task File
- If creating: Use `write_to_file` to write the complete task document.
- If updating: Use `replace_file_content` or `write_to_file` (with overwrite) to update the existing task file with the refreshed cumulative state.

### Step 4: Report to User
Provide a concise confirmation linking to the task file:
- Target file path (e.g., [`.agents/tasks/10_...md`](file:///...)).
- Brief bulleted summary of what was recorded.
