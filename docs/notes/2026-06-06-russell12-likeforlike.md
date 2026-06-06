# Russell-12 like-for-like closure (#135) — the decisive diagnostic

Date: 2026-06-06. Follow-up to the #125 Part-2 campaign
(`docs/notes/2026-06-06-validation-campaign-125.md`). Golden discipline absolute:
the EXPECTED side of every match check and every seed is the row's SOURCED
anchor (Russell 2004 / McConaghy values, segment ToFs, descriptor arc ToFs).
Nothing our own code computes is ever used as an EXPECTED; nothing loosened.

## The question

The #125 Part-2 sweep drove the N-arc ballistic corrector on
`Ephemeris('astropy')` (real DE440) and compared closures to the rows'
CIRCULAR-COPLANAR sourced anchors. Cross-fidelity confounded the model mismatch
with the basin problem (0 CLOSE-AND-MATCH, 11 off-anchor, vector mode 0).

But every Russell/McConaghy row is, by construction, a solution OF the
circular-coplanar model. So the test that removes the confound is: run the SAME
genome derivation, tolerances and anchors on `Ephemeris('circular')` —
coplanar-vs-coplanar, like-for-like, against known-solvable instances.

Driver: `scripts/campaign_russell12.py --model circular` (the #125 runner,
extended to thread the ephemeris token through the worker and add the
seed-at-truth probe; the default `--model astropy` path is byte-unchanged so
#125 still reproduces). 256-epoch grid over one target period, both residual
modes, 16 workers. Tolerances unchanged: V∞ 0.5 km/s, transits 5 d,
period 0.05 yr, corrector floor 0.1 km/s.

## Per-row table (verbatim, model=circular)

| id | outcome | seq | revs | closed mag/vec | E∞ src/ach | M∞ src/ach | transits src/ach |
|---|---|---|---|---|---|---|---|
| mcconaghy-2006-em-k2  | CLOSE-OFF-ANCHOR | E-M-E-E   | 0,0,2   | 207/0 | 4.70 / 8.97  | 5.00 / 9.84  | 153,153 / 168.9,168.9 |
| russell-ch4-4.991gG2  | CLOSE-OFF-ANCHOR | E-M-E-E   | 0,0,2   | 201/0 | 4.99 / 8.97  | 5.10 / 9.84  | 150,150 / 168.9,168.9 |
| russell-ch4-8.049gGf2 | CLOSE-OFF-ANCHOR | E-M-E-E   | 0,0,1   | 127/0 | 8.05 / 11.87 | 10.02 / 16.20| 93,93 / 63.7,115.3 |
| russell-ch4-9.353Gg2  | CLOSE-OFF-ANCHOR | E-M-E-E   | 0,0,2   | 120/0 | 9.35 / 8.97 OK | 10.52 / 9.84 | 85,85 / 168.9,168.9 |
| russell-ch4-3.64gGg3  | NO-CLOSE         | E-M-E-E-E | 0,0,2,1 | 0/0   | 3.64 / —     | 4.59 / —     | — |
| russell-ch4-3.78Gg3   | CLOSE-OFF-ANCHOR | E-M-E-E-E | 0,0,3,1 | 231/0 | 3.78 / 23.85 | 4.63 / 14.79 | 171,171 / 209.0,156.3 |
| russell-ch4-5.30gGf3  | CLOSE-OFF-ANCHOR | E-M-E-E-E | 0,0,1,3 | 211/0 | 5.30 / 23.85 | 9.17 / 14.79 | 118,118 / 209.0,156.3 |
| russell-ch4-9.94Gg3   | CLOSE-OFF-ANCHOR | E-M-E-E-E | 0,0,4,1 | 19/0  | 9.94 / 28.05 | 10.76 / 18.40| 82,82 / 234.8,234.8 |
| russell-ch4-3.66gfF3  | NO-CLOSE         | E-M-E-E-E | 0,0,1,3 | 0/0   | 3.66 / —     | 4.66 / —     | — |
| russell-ch4-5.30ggF3  | CLOSE-OFF-ANCHOR | E-M-E-E-E | 0,0,1,3 | 221/0 | 5.30 / 23.85 | 5.44 / 14.79 | 143,207 / 59.0,306.3 |
| russell-ch4-5.75ggF3  | NO-CLOSE         | E-M-E-E-E | 0,0,2,2 | 0/0   | 5.75 / —     | 9.36 / —     | — |
| russell-ch4-6.44Gg3   | CLOSE-OFF-ANCHOR | E-M-E-E-E | 0,0,3,1 | 91/0  | 6.44 / 9.55  | 3.74 / 5.35  | 262,262 / 243.5,243.5 |

(The "ach" column is OUR computed V∞ at the closed geometry — off-anchor
evidence, never asserted as an EXPECTED. Period pins to the sourced value on
every closing row by construction of the slack leg.)

### Counts (model=circular)

- **CLOSE-AND-MATCH: 0**
- **CLOSE-OFF-ANCHOR: 9**
- **NO-CLOSE: 3** (`russell-ch4-3.64gGg3`, `russell-ch4-3.66gfF3`, `russell-ch4-5.75ggF3`)

## The headline

**Like-for-like did NOT change the verdict: still 0 CLOSE-AND-MATCH.** The model
mismatch was not the confound. Running coplanar-vs-coplanar against rows that are
provably solutions OF the coplanar model, the descriptor-seeded corrector still
lands in the high-V∞ degenerate basin (E∞ 8.97–28 km/s vs sourced 3.6–9.9;
transits drift off the symmetric sourced ToFs) and vector mode still closes 0
(the degenerate basin is bend-infeasible). Three rows that closed off-anchor on
astropy are NO-CLOSE on circular — the circular Lambert geometry is *less*
forgiving at the seeded topology, not more.

## The seed-at-truth probe (the diagnosis)

For each row, the probe fixes the free-leg ToFs at the row's OWN sourced geometry
(transit ToFs from `trajectory.segments`, E-E loop ToFs from
`free_return_arcs[]`), scans t0 over one period to find the phase that minimises
the residual evaluated AT that truth geometry, then seeds the corrector EXACTLY
there and asks: does it STAY (residual→0 at truth, ToFs unmoved) or WALK AWAY?

| id | best-phase truth residual km/s | truth below 0.1 floor? | solved residual | ToF drift d | verdict |
|---|---|---|---|---|---|
| mcconaghy-2006-em-k2  | 3.43   | NO | 0.000 | 15.9  | WALKED-AWAY |
| russell-ch4-4.991gG2  | 4.00   | NO | 0.000 | 18.9  | WALKED-AWAY |
| russell-ch4-8.049gGf2 | 3.18   | NO | 0.000 | 12.9  | WALKED-AWAY |
| russell-ch4-9.353Gg2  | 10.42  | NO | 8.301 | 19.7  | WALKED-AWAY |
| russell-ch4-3.64gGg3  | 1000.0 | NO (Lambert pathology) | inf | 0.0 | WALKED-AWAY |
| russell-ch4-3.78Gg3   | 18.04  | NO | 0.000 | 38.0  | WALKED-AWAY |
| russell-ch4-5.30gGf3  | 10.72  | NO | 0.000 | 29.0  | WALKED-AWAY |
| russell-ch4-9.94Gg3   | 17.16  | NO | 5.634 | 107.9 | WALKED-AWAY |
| russell-ch4-3.66gfF3  | 1000.0 | NO (Lambert pathology) | inf | 0.0 | WALKED-AWAY |
| russell-ch4-5.30ggF3  | 37.55  | NO | 0.000 | 156.9 | WALKED-AWAY |
| russell-ch4-5.75ggF3  | 1000.0 | NO (Lambert pathology) | inf | 0.0 | WALKED-AWAY |
| russell-ch4-6.44Gg3   | 17.16  | NO | 0.000 | 32.1  | WALKED-AWAY |

**All 12 WALK AWAY — and the reason is one level deeper than seeding.** On NO row
is the residual at the sourced geometry anywhere near the 0.1 km/s closure floor:
the *smallest* best-phase truth residual is 3.18 km/s (8.049gGf2), the typical
value is 10–37 km/s, and three rows (the `gfF`/`ggF`/`gGg` topologies) hit the
Lambert-pathology penalty (1000) at truth — no t0 phase makes the sourced ToFs
even Lambert-solvable in our genome.

In other words: **the sourced geometry is not a closure point of our genome.**
The corrector does not abandon a good seed for a degenerate one (that would be a
seeding/solver-attraction story); it correctly walks away from a point that
genuinely does not satisfy our V∞-continuity residual, toward one that does.

## A/B verdict

This is **neither clean Outcome A nor the simple Outcome B** of the brief.

- Not A: no like-for-like match exists; there is no multi-arc validation evidence
  to report, and the §14-V1 like-for-like question is moot for these rows (there
  is no reproduction to grade).
- Not the simple B ("family selection is our seeding/corrector deficiency, the
  basin exists but we don't seat it"): the seed-at-truth probe shows the sourced
  basin **does not exist in our genome as posed**. Seeding at truth cannot help
  because truth is not a residual-zero point of our model.

The deficiency is in the **genome mapping** (descriptor → corrector free
variables), not in seeding and not in the least-squares solver. Candidate
culprits, in order of suspicion:

1. **Transit-leg modelling.** The two Mars transfer legs are seeded as the
   sourced single-rev E→M and M→E ToFs, but the corrector then treats them as
   free single-ellipse Lambert arcs between the circular planet positions. The
   sourced 153 d symmetric transfer is a *specific* low-V∞ ellipse; the Lambert
   arc our genome solves at the same ToF, between the circular E and M positions
   at the corrector's chosen phase, is a different, higher-V∞ ellipse. The 3–4
   km/s truth residual on the k=2 rows is exactly this transfer-geometry gap.
2. **E-E loop rev/branch assignment.** The `per_leg_revs` bracketing
   (`round(seed/yr) − 1`) and the fixed `"low"` branch may not be the rev/branch
   the sourced loop actually flies; a wrong rev count makes the loop Lambert
   either pathological (the three 1000-residual NO-CLOSE rows) or a different
   ellipse.
3. **Slack-leg elimination.** Pinning the period by reconstructing the longest
   E-E loop as `period − Σfree` couples the period constraint into a single loop;
   if that loop's true ToF is not the residual, the reconstructed ToF forces an
   off-geometry Lambert.

The right next step is NOT a denser scan or a better seed ladder — those cannot
reach a basin that is not there. It is to re-derive the genome so that the
sourced geometry IS a residual-zero point: most likely solving the transfer legs
as constrained (a, e) arcs matching the sourced V∞ rather than free Lambert arcs
at fixed ToF, and deriving each E-E loop's rev/branch from its sourced ToF +
resonance rather than bracketing. This is squarely an M-ED genome-modelling item,
not an M-3D / real-eccentricity one (the like-for-like result rules out fidelity
as the cause).

## Golden / writeback discipline

No catalogue field is touched. No `validation_level` writeback (there is no match
to promote, and no reproduction to grade against §14 V1). The decisive artefact
is the probe table above; the campaign JSON is `/tmp/russell12_circular.json`
(regenerable with the command below). A pinning test
(`tests/test_russell12_likeforlike_probe.py`) freezes the headline invariants:
0 CLOSE-AND-MATCH on circular, and the seed-at-truth residual floor being
unreachable on the representative symmetric row.

Reproduce::

    uv run python scripts/campaign_russell12.py --model circular \
        --epochs 256 --workers 16 --probe-at-truth --phase-epochs 256 \
        --out /tmp/russell12_circular.json
