# ER3BP discovery negatives — literature erratum (2026-06-25)

**Trigger:** acquisition + digest of five ER3BP/3D periodic-orbit references
(2026-06-25) that we lacked when #432/#435/#436 were closed. Three of them
directly bear on the ER3BP e>0 negatives. This note records how each negative
must be re-scoped. **No campaign result is overturned silently; each is narrowed
or reclassified with its source.**

## The three ER3BP references

1. **Antoniadou & Libert 2018**, "Origin and continuation of 3/2, 5/2, 3/1, 4/1,
   5/1 resonant periodic orbits in the circular and elliptic RTBP," CMDA, DOI
   10.1007/s10569-018-9834-8, arXiv:1805.00288. (digest:
   `docs/notes/2026-06-25-digest-antoniadou-libert-2018-circ-ellip-continuation.md`)
2. **Peng, Bai & Xu 2017**, "Continuation of periodic orbits in the Sun-Mercury
   elliptic restricted three-body problem," Commun. Nonlinear Sci. Numer.
   Simulat. 47:1-15, DOI 10.1016/j.cnsns.2016.11.005. (digest:
   `2026-06-25-digest-peng-2017-sun-mercury-ERTBP.md`)
3. **Martínez-Cacho, Gil Calvo, Bombardelli & Baresi 2025**, "Planar retrograde
   periodic orbits in the elliptic restricted three-body problem," Acta
   Astronautica 229:430-465, DOI 10.1016/j.actaastro.2025.01.006. (digest:
   `2026-06-25-digest-planar-retrograde-ERTBP-2025.md`)

## The core finding — our e>0 negatives are (partly) method-blind-spot

**Antoniadou & Libert 2018 prove ISOLATED high-eccentricity elliptic families
exist** — "new isolated (not continued from bifurcation points) families... at
high eccentricities," one per resonance (3/2, 5/2, 3/1, 4/1, 5/1), predominantly
stable, in fixed symmetric configs. **By construction these have no circular
limit, so any continuation-FROM-circular method is structurally blind to them.**
#432 (continue known cyclers into e) and #436 (reverse-continue to e=0) are
exactly such methods. Therefore:

> The #432/#436 "no novel e>0-only family" results are **method-blind-spot
> negatives, not non-existence proofs.** The published ER3BP record contains a
> family class our methods cannot reach.

**Caveat that keeps this honest:** Antoniadou & Libert's isolated families are
*eccentric mean-motion-resonant* families at µ≈1e-3 (Sun-giant), **not Earth-Moon
and not co-orbital cyclers.** So this does NOT prove Earth-Moon cycler isolated
families exist — it proves the *phenomenon* exists in the ER3BP, which removes our
license to assert "none exist." Whether the Earth-Moon cycler class has isolated
e>0 members is now genuinely **open**, and is the subject of task #440 (direct
high-e seeding).

## Per-campaign re-scoping

### #435 — split verdict
- **DRO/retrograde leg: CORROBORATED.** Martínez-Cacho 2025 (QSO = DRO = Hénon
  *f* family) shows every retrograde co-orbital family is circular-rooted
  (Hill→CR3BP→ER3BP); no isolated elliptic-only steady DRO. Our #435 DRO negative
  stands, now with a published same-model backing.
- **Lyapunov leg: NARROW to its actual scope.** Peng 2017, in the *same*
  Sun-Mercury system at the *same* e=0.206, continues *multi-revolution resonant
  Halo* families and finds non-monotonic-in-e characteristic curves with folds,
  disconnected branches (A/B/C), an isolated loop, and complex-instability
  transitions — and parts that **cannot reach** the real Mercury eccentricity.
  #435 used *low-revolution Lyapunov/DRO* seeds and saw none of this. So #435's
  "survive smoothly, no novel family" is a property of the **low-rev seed
  choice**, not of the system. Re-scope #435's claim to "low-rev Lyapunov/DRO
  seeds under naive continuation," and cite Peng as the same-system multi-rev
  counterexample. (No isolated elliptic-only family in Peng either — all his
  families attach to e=0 — so #432/#436's isolated-family gap is unchanged.)

### #432 / #436 — reclassify as method-blind-spot (see core finding). #436's
own verdict already flagged its discriminator unreliable (classification flipped
on step count); we now also know its *targets* (isolated families) are
unreachable by continuation at all. Both effects point to the same fix.

### #437 (fold-aware pseudo-arclength) — VALIDATED, with goldens
Both new papers independently confirm the #436 failure mode and #437's fix:
- Peng: "we tested to directly continue... the routine failed" at turning points;
  fix = pseudo-arclength/tangential continuation + multiple shooting. **Tables
  2 & 3** give eight exact Sun-Mercury M5N2 ME-Halo ICs + monodromy eigenvalues
  — a same-µ same-e sourced golden (ER3BP has no Jacobi integral; match the
  eigenvalue triple).
- Martínez-Cacho: explicit **fold bifurcation at e=0.0324** (Sun-Mars 2:3 swing
  QSO), Appendix B 15-digit QSO ICs (CR3BP + ER3BP), with the Sun-Mars 2:3 ER3BP
  entry left **blank** ("couldn't close") — the fold-failure regression golden.

## Disposition / follow-ups
- **#437** now has two sourced same-model goldens (Peng Tables 2/3; Martínez-Cacho
  App. B + the e=0.0324 fold). Build the fold-aware continuator against these.
- **#440** (new): direct high-e seeding of isolated families in Antoniadou &
  Libert's symmetric configs (their µ=0.001 same-model golden: gap a1 =
  0.763143 / 0.542884 / 0.39685, bifurcation periods T0 = 4π, 4π/3, 2π, 2π/3, 2π).
- **Registry:** the `er3bp-*` entries in `data/empty_regions.jsonl` carry a
  pointer to this erratum (re-scoped, not deleted — they remain valid statements
  about what each *method* found).
- **New acquisition targets surfaced by the digests:** Breakwell & Brown 1979
  (Celest. Mech. 20, 389) — Earth-Moon halo ICs (Howell's seed source; the
  Earth-Moon 3D coordinate catalogue we still lack); Ollé & Pacha 1999
  (bifurcation from limiting restricted problems, Peng's cited mechanism);
  Peng & Xu 2015 CMDA (prior ME-Halo stability).
