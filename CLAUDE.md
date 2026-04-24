# CLAUDE.md

## Commands

### /handoff

When the user says "/handoff" or "hand off":

Write `HANDOFF.md` in the project root with everything a fresh agent needs to finish the current task. Include:

1. **Goal** — What we're trying to accomplish (the original ask)
2. **What's done** — Changes made so far, files touched, decisions locked in
3. **What worked** — Approaches/tools/patterns that succeeded
4. **What didn't work** — Failed approaches and why they failed (so the next agent doesn't repeat them)
5. **What's left** — Remaining steps to complete the task
6. **Key context** — Non-obvious gotchas, constraints, or domain knowledge discovered along the way

The file must be **self-contained**. A new agent loading only HANDOFF.md should have enough context to continue and finish without asking questions about prior attempts.
