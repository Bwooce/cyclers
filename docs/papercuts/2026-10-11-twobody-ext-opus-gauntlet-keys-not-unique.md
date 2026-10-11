# The two-working-body gauntlet JSONs can hold two candidates with the same "key"

- Date: 2026-10-11
- Agent: twobody-ext-opus
- Seen before: no (it bit #1024's first continuation run; #1025's realeph lookup for ev members used key + V_inf and was not affected)

What happened: a `*_gauntlet.json` candidate's `key` is the cycle key of its best member. Two merged physical groups (different V_inf/turn sequences) can share that key, because one structure can have several zeros. In #998's Pluto gauntlets, keys repeat (ps 183 candidates, 181 distinct keys; pk 62, 56). #1024's first run keyed its records by `key`, which silently lost members and mis-matched the realeph lookup.
Workaround: #1024 keys records by gauntlet index, and matches by key plus ideal worst ratio.
Suggested fix: give every gauntlet candidate a stable id (cell + index, or key + a hash of x_days), and use it in every downstream lookup. Grep the scripts that build dicts keyed by `c["key"]` from a gauntlet file.
