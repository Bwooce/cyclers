# Window-4 apply scripts kept in the session scratchpad were lost in a session restart

- Date: 2026-10-12
- Agent: team-lead (cyclers)
- Seen before: no

What happened: two agents built the window-4 catalogue apply scripts (gc4-19 row, six Earth-Venus V2 promotions) in their scratchpad folders and reported path and md5 to the lead. The session process ended before the window opened; the scratchpad was wiped with it, and both scripts were gone. Their inputs were committed, so a fresh agent is rebuilding them (about an hour of work repeated).
Workaround: respawned an agent to rebuild both scripts from the committed drafts and judge files.
Suggested fix: any artifact that another agent or the lead must consume later (apply scripts, patch files, review drafts) lives in the repository (committed, or at least untracked in a repo path such as scripts/window4/), never only in the scratchpad. Rule added to docs/team/coordination.md.
