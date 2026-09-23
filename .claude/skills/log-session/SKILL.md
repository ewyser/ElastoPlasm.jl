---
name: log-session
description: Reconcile ElastoPlasm.jl's doc set (CLAUDE.md, .claude/docs/, .claude/bug/known/, .claude/bug/fixed/) with what actually happened in this conversation — fixed bugs, newly found bugs, design decisions, follow-up work. Invoke at the end of a work session, or whenever asked to "write this up" / "update the docs" / "log what we did".
---

# Log session

This repo's working memory is split across `CLAUDE.md` (overview), `.claude/docs/*.md`
(architecture/operations/persistence/gotchas/conventions/planned-improvements), and
`.claude/bug/known/`+`.claude/bug/fixed/` (one file per bug). This skill's job is to fold this
session's real work into that structure correctly — not to write a session transcript
or a new standalone report file.

## Step 1 — Establish what actually happened

Re-read the conversation and classify each substantive thing that happened into
exactly one of these buckets. Do not log a bucket that didn't happen — an empty
session on some axis means no edit to that file, not a filler sentence.

- **A bug was fixed.** A `.claude/bug/known/*.md` file's issue was resolved, or a new bug was
  found *and* fixed in the same session.
- **A new bug was found but not fixed.** Needs a new `.claude/bug/known/*.md` file.
- **A design/architecture decision was made or changed** (a new type, a renamed
  convention, a new config knob, a dispatch pattern). Belongs in
  `.claude/docs/architecture.md`, `.claude/docs/conventions.md`, or
  `.claude/docs/planned-improvements.md`, whichever already covers that topic.
- **Follow-up work was identified but deliberately not done.** Belongs in
  `.claude/docs/planned-improvements.md` as a new bullet, not a bug file — this
  distinction already matters in this repo (see that file's own framing: "not bugs,
  just known follow-up work").
- **Nothing doc-worthy happened** (pure exploration, a question answered from
  existing docs, a false lead that was ruled out with no lasting conclusion). Say so
  and stop — don't force a write.

## Step 2 — Match this session's writing conventions

Before writing anything, skim 2-3 existing files in the bucket you're writing to
(e.g. two files under `.claude/bug/fixed/` if closing a bug) and match their shape exactly:

- `.claude/bug/*/*.md`: one bug per file, kebab-case filename from the lead symptom, a
  `# <symptom>` heading, a `**Status: fixed.**`/`**Status: open.**` line, then the
  narrative — root cause, fix, and any methodology lesson or follow-on finding,
  written the same terse, evidence-citing style already there (concrete numbers,
  concrete file/line references, not vague summary).
- `.claude/docs/*.md`: no per-file frontmatter, just a `##`-level section appended or
  edited in place inside the existing file — do not create a new doc file for a topic
  one of the six existing files already owns.
- Cross-references: use the same relative-path style already in use
  (`` `.claude/bug/known/foo.md` ``, `` `.claude/docs/architecture.md` ``) — check both
  directions: a new bug file that supersedes/relates to an existing one should link
  it, and the existing one should be updated to point forward if it's now stale or
  fixed.
- If a bug moves from open to fixed: move the file from `.claude/bug/known/` to `.claude/bug/fixed/`
  (git mv, don't leave a duplicate), update its status line, and grep the rest of the
  tree for anything that referenced it by its old path.

## Step 3 — Apply the edits

Write/move the files. Keep changes scoped to documentation — this skill does not
touch source code. If source code also needs a change as a result of this session's
findings, say so explicitly and ask before making it (or hand off to the user), don't
fold a source fix into a "log session" pass silently.

## Step 4 — Report back

End with a short summary: which files were added/edited/moved, and one line per file
on what changed — not a repeat of the full file contents.
