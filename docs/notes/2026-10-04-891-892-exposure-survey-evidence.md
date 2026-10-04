# Exposure of past results to the #891 / #892 model defects: file and line evidence

Appendix to `2026-10-04-891-892-dependent-code-triage.md` (section 5), written by the triage agent's
read-only survey on 2026-10-04 and preserved here by the coordinator so the evidence outlives the
session. Line numbers refer to the tree at the commit named below.

Read-only survey, HEAD at survey time 69184502. Nothing was run.
JL = line number in data/empty_regions.jsonl.

## Importer map (what loads the defective modules)

Direct importers of core.bcr4bp (grep of `src`):
- genome/bcr4bp_genome.py:57, bcr4bp_continuation.py:47, bcr4bp_systems.py:55, bcr4bp_torus.py:18, bct_transfer.py:43
- core/wsb.py:61 (calls bcr4bp.propagate_bcr4bp at wsb.py:240), search/cislunar_bct_search.py:36 (+ wsb:37, bct_transfer:38),
  search/sun_forced_periodic_884.py:78, search/bvp_integral.py:53
- data/validation/v0_bcr4bp.py:56, v1_bcr4bp.py:72, v2_bcr4bp.py:60, v3_bcr4bp.py:66
Direct importers of core.qbcp: genome/qbcp_torus.py:19, search/variational_qbcp_torus.py:177,
  search/variational_qbcp_arc.py:197, search/variational_periodic_orbit_qbcp.py:106.
Imported-but-only-helpers (not exposed through physics): search/variational_ccr4bp_torus.py:112 and
  search/variational_crnbp_torus.py:90 import `_basis_matrices` / `_k2_first_harmonic_cols` from variational_qbcp_torus; the file's own comment
  (variational_ccr4bp_torus.py:108-111) says these depend only on mode counts, not physics. `grep 'qbcp\.'` in both files: no hits.
  So the ccr4bp/crnbp torus scripts (screen_694/695/696/701/703/716/882, run_704/705/729, verify_721/724) are NOT exposed (consistent with OUTSTANDING.md:~968).
data/validation/__init__.py imports v0..v3_bcr4bp (lines 12/24/36/66), so merely importing the package loads the module; only callers of run_v*_bcr4bp use it.

## (i) data/empty_regions.jsonl entries mentioning bcr4bp / bicircular (19 lines)

Verified count: 19 lines (JL 25, 36, 39, 51-58, 59, 60, 61, 62, 63, 81, 83, 94).
Check method for every row: list of scripts under scripts/ importing bcr4bp/qbcp/bct/wsb/sun_forced (script_imports.txt in this directory) was compared with the
producer named in the entry's `run` / `source_anchors` fields.

| JL | region_id | producer (from entry) | verdict | evidence |
|---|---|---|---|---|
| 25 | er3bp-discovery-em-broucke-koblick-2026-06-24 | scripts/run_432_er3bp_discovery.py (JL25 run.script; file exists) | PROSE | Only mention is the scoping-note filename `2026-06-16-frontier-scoping-er3bp-bcr4bp-3d-qp-epoch.md` in source_anchors (JL25). Method is ER3BP (method_capability.genome "#293 ER3BP genome", JL25). run_432 is not in the bcr4bp/qbcp importer list; genome/er3bp_*.py do not import bcr4bp. |
| 36 | isolated-er3bp-cycler-novelty-gap-analysis-2026-06-25 | none: run = {task 442, "pre-build novelty-gap adjudication"} (JL36) | PROSE | Mention is the same scoping-note filename in a citation list (JL36). No script, no integration. |
| 39 | cislunar-bct-wsb-quasicycler-2026-06-26 | #378; run.git_sha 6cd486d, run.note docs/notes/2026-06-26-378-cislunar-bct-verdict.md (JL39). Writer scripts/_apply_378_empty_region.py (exists; only appends the report, imports no model). Sweep code search/cislunar_bct_search.py; apoapsis spike scripts/spike_378_bct_apoapsis.py | INTEGRATED | JL39 method_capability.genome: "backward BCT constructor (genome/bct_transfer.construct_bct_backward) on the incoherent BCR4BP (core/bcr4bp...)"; centre: "Earth-Moon BCR4BP with the Andreu Sun perturbation". Code: search/cislunar_bct_search.py:36 imports core.bcr4bp, :37 core.wsb, :38 genome.bct_transfer; core/wsb.py:61,240 propagates bcr4bp; bct_transfer.py:43; scripts/spike_378_bct_apoapsis.py:35 imports core.bcr4bp. Note the verdict note (docs/notes/2026-06-26-378-cislunar-bct-verdict.md:15,31) names cislunar_bct_search as the 44-point grid sweep. Could not tell which script under scripts/ launched the 44-point grid (no script in scripts/ imports cislunar_bct_search); the module is the producer. |
| 51 | floquet-branch-c32-b0-em-v4-2026-06-19 | #389, git_sha 27a3959 (JL51 run) | PROSE | bcr4bp appears only in a #425 `reverification` text: "CR3BP/BCR4BP have no epoch & no flyby" (JL51). Method is real-ephemeris / ballistic closure; no script importing bcr4bp is a floquet-branch or V4 closer. |
| 52 | floquet-branch-c32-c-3.1774-em-v4-2026-06-19 | #392, 9286782 | PROSE | same reverification sentence (JL52) |
| 53 | floquet-branch-c11a-b0-em-v4-2026-06-19 | #393, 9286782 | PROSE | same (JL53) |
| 54 | floquet-branches-em-cycler-nodes-v4-aggregate-2026-06-19 | #392, 9286782 | PROSE | same (JL54) |
| 55 | closer-sweep-v1-russell-ocampo-3.1.2+1-2026-06-17 | #365, fd95838 | PROSE | same sentence (JL55); heliocentric Russell-Ocampo closer sweep, not an EM model |
| 56 | closer-sweep-v1-russell-ocampo-4.3.1-5-2026-06-17 | #365 | PROSE | same (JL56) |
| 57 | closer-sweep-v1-russell-ocampo-4.5.2-2-2026-06-17 | #365 | PROSE | same (JL57) |
| 58 | closer-sweep-v1-mcconaghy-2006-em-k2-s1l1-2026-06-17 | #365 | PROSE | same (JL58) |
| 59 | cross-system-se-em-l2-patched-cr3bp-2026-06-20 | #405, 452b28d; scripts/run_405_cross_system_search.py, analyze_405_theta_closure.py (exist) | PROSE | BCR4BP appears as an untested next venue: "(2) BCR4BP with an SE-scale seed (#412 re-scope)" (JL59) plus the #425 reverification sentence. Method_capability tags: cr3bp, patched-cr3bp (JL59). genome/cross_system_cycle.py has no bcr4bp/qbcp reference (grep) and neither run_405 nor analyze_405 is in the importer list. |
| 60 | bcr4bp-phase-b-em-libration-seed-reach-spike-2026-06-20 | #412, git_sha 3f51739; scripts/spike_412_bcr4bp_reach.py (exists); note docs/notes/2026-06-20-412-bcr4bp-reach-spike.md | INTEGRATED | JL60 method_capability: "CR3BP EM-L1 Lyapunov orbit cast as a BCR4BP-at-mu_sun=0 seed (correct_bcr4bp_periodic)", corrector "continue_bcr4bp_family_in_musun". scripts/spike_412_bcr4bp_reach.py:21 imports core.bcr4bp, :23 genome.bcr4bp_continuation, :24 genome.bcr4bp_genome, :33 genome.bcr4bp_systems. The Sun enters the continuation, so the reach numbers (0.9-1.16 LD) were measured in the wrong-sense model. |
| 61 | cross-system-se-em-3d-patched-cr3bp-2026-07-01 | #515, c65f858; scripts/run_515_cross_system_3d_search.py | PROSE | "RESWEEP CONDITION: BCR4BP 3D halo continuation ..." is future work (JL61); tags patched-cr3bp (JL61); run_515 not in importer list |
| 62 | cross-system-se-em-3d-asymmetric-patched-cr3bp-2026-07-01 | #517; scripts/run_517_asymmetric_3d_search.py | PROSE | same pattern (JL62) |
| 63 | cross-system-se-em-3d-multirev-patched-cr3bp-2026-07-01 | #516; scripts/run_516_multirev_3d_search.py | PROSE | same pattern (JL63) |
| 81 | mars-phobos-deimos-symmetric-closure-609-2026-07-16 | #609; scripts/_apply_609_mars_phobos_deimos_empty_region.py + data/enumerate_609_*.jsonl (JL81 run.data_files) | PROSE | Mention is "no published Sun-Mars-Phobos BCR4BP ... cycler exists", a literature-gap remark about literature_check.py's anchor (JL81). Neither the _apply_609 script nor the enumeration is in the importer list. (src/cyclerfinder/search/literature_check.py does contain the string bcr4bp; it is a signature/anchor table, not a propagation - not read in detail, so "could not tell" for any literature_check internals, but this entry's method is a symmetric-closure enumeration.) |
| 83 | cross-system-se-em-l1l1-patched-cr3bp-2026-07-17 | #622; scripts/run_622_em_l1_se_l1_search.py (JL83 run.script; exists) | PROSE | "BCR4BP with an SE-scale seed (#412 re-scope)" is the untried alternative (JL83). run_622 not in importer list. |
| 94 | sunmars-bct-wsb-quasicycler-2026-07-22 | #681; scripts/reproduce_681_topputo_belbruno.py, scripts/search_681_sunmars_wsb_chain.py (named in source_anchors, JL94) | PROSE | Method_capability.genome (JL94): new core/sunmars_wsb.py "physical planar Sun-Mars ER3BP integrator in a heliocentric inertial frame", with the sentence "Generalizes the #378 cislunar W-set predicate (which was Earth-Moon-BCR4BP-specific)". grep: core/sunmars_wsb.py has no bcr4bp/qbcp reference; search_681 imports core.sunmars_wsb (line 59), not core.bcr4bp/wsb. Unaffected, though note the concept lineage (#378). |

Count: INTEGRATED = 2 (JL39, JL60); PROSE = 17; UNKNOWN = 0.
OUTSTANDING.md:967-969 says "19 negative-result stamps ... mention the bicircular model". That is the grep count; only 2 of the 19 rest on an
integration of the model. Rows JL51-58 and 59 additionally carry the generic #425 reverification sentence which concerns "CR3BP/BCR4BP have no epoch & no flyby" as a
structural argument (a statement about model scope, unaffected by the Sun sense). Recommended stamping: METHOD-INVALID for JL39 and JL60 only.

Caveat on JL39: the #378 negative was about a W-set / ballistic-capture quasi-cycler in the incoherent BCR4BP; JL39 itself lists as re-open keys
"(1) forward-from-LEO on-W capture did not converge in the incoherent BCR4BP" and "(2) the coherent QBCP". Both are model-bound.

## (ii) Task conclusions and dependence on the defective model

Source: data/OUTSTANDING.md and docs/notes/. "Yes" = the conclusion's evidence was computed by integrating core/bcr4bp.py or core/qbcp.py as defective.

| Task | What it concluded | Depends on defective model? | Evidence |
|---|---|---|---|
| #292 | BCR4BP Phase 1 substrate landed: EOM, STM, propagator, periodic-orbit corrector; 17 tests; CR3BP limit (mu_sun=0) holds to float precision; "POL1 seed closes to a nearby BCR4BP L1 dynamical substitute"; no catalogue writeback, no novelty claim. | Partly. It IS the defective module. The CR3BP-limit, STM-vs-finite-difference and constants checks are sense-independent (the Sun term vanishes or the check is self-referential), so those stand. The Sun-on closure claim (POL1 seed substitute) and the Sun-position kinematics test are in the non-physical model. | docs/notes/2026-06-16-292-bcr4bp-phase1.md:7-12,31-37 (status, tests); OUTSTANDING.md:967 lists #292 as exposed; OUTSTANDING.md:~960 "no test against a published bicircular ORBIT ... every later check compared the model with itself" |
| #303 | Phase 2: mu_sun natural-parameter continuation of the CR3BP L1 Lyapunov (C=3.1294) to Andreu mu_sun=328900.54; 50/50 steps converged; IC barely moves (1e-4 in x, 1e-2 in vy, T 2.946->2.950 TU); stays hyperbolic_pair; catalogue probe 0 candidates / 11 no-match. | Yes. Every Sun-on number (final x0, vy0, T, stability tag, the 11 no-match records) comes from the wrong-sense model. Outputs data/bcr4bp_l1_family_303.jsonl (51 lines), data/bcr4bp_validation_bridges_303.jsonl. | docs/notes/2026-06-16-303-bcr4bp-phase2-musun-continuation.md:25-60; scripts/run_303_bcr4bp_l1_continuation.py:37,56,270; OUTSTANDING.md:967 |
| #304 | Phase 3: halo corrector mask + Howell EM-L1 southern halo continued in mu_sun to Andreu value; 50/50 converged; z0 -0.0529 -> -0.0566, T 2.760 -> 2.777 TU; hyperbolic_pair; probe 0 candidates / 11 no-match. | Yes. Same reasoning; outputs data/bcr4bp_halo_family_304.jsonl, data/bcr4bp_halo_validation_bridges_304.jsonl. | docs/notes/2026-06-16-304-bcr4bp-halo-phase3.md:27-50; scripts/run_304_bcr4bp_halo_continuation.py:39,53,303; OUTSTANDING.md:967 |
| #334 | Phase 4: swept the L1 Lyapunov mu_sun-continuation over 8 Sun-primary-secondary systems; fitted scaling dx0_target/dx0_SEM ~ (mu_sun ratio)(a_sun ratio)^k, k = 2.89 +- 0.27 inside the #326 band [2,3]; flagged Pluto-Charon outlier and Sun-Saturn-Titan high ratio; added 4 KNOWN_CORPUS anchors; no writeback, no novelty claim. | Partly. The fitted k and the outlier flags are from Sun-on continuations in the wrong-sense model (yes). The constants table (derived omega_sun, mu_sun, a_sun; note's section 1) and the 4 KNOWN_CORPUS literature anchors do not depend on the Sun's sense (no). Could not tell whether k survives; the displacement magnitude scales with the Sun's tidal strength, which the sense error may not change, but this is untested. | docs/notes/2026-06-17-334-bcr4bp-system-swap.md:6-22 (TL;DR), 24-60 (constants); scripts/scan_334_bcr4bp_system_swap.py:80-92,103,666; OUTSTANDING.md:967 |
| #412 | Phase-B reach spike: EM-L1 Lyapunov seeds (C=3.1294/3.05/2.95) cast into the BCR4BP and mu_sun-continued keep Earth-relative reach ~0.9-1.16 LD vs the 3.90 LD SE target; family is the wrong vehicle for the cross-system cycle; scoped negative; re-scope needs an SE-scale BCR4BP seed. | Yes. The reach numbers and the "family does not survive to full mu_sun" statements are outputs of the wrong-sense continuation (JL60). The structural argument (an EM-bounded orbit stays EM-scale) may still hold; could not tell. This is one of two INTEGRATED stamps. | data/empty_regions.jsonl:60 (result.reach_by_C); scripts/spike_412_bcr4bp_reach.py:21-33; docs/notes/2026-06-20-412-bcr4bp-reach-spike.md:1-18; OUTSTANDING.md:14171 |
| #533 | Built the "genuine coherent QBCP" model (8 Fourier alpha series, Gimeno-Jorba 2018 Table 4), EOM + STM in canonical variables; "verified against circular BCR4BP limits and finite differences" (resolved 2026-07-08). | Yes. The model itself had wrong parities and a misplaced alpha_6 (#892, OUTSTANDING.md:~872-915); its own verification (BCR4BP limit + finite differences) could not catch either (593 note says FD verifies a formula against its own derivative). SE-L2 torus numbers quoted as the clean baseline are not reproducible. | OUTSTANDING.md:7391; docs/notes/2026-07-14-593-qbcp-alpha6-impact-scoping.md:63,78 ; OUTSTANDING.md:888,912 (#892 lists #533 as exposed); scripts/run_533_qbcp_connection.py:20,22 |
| #538 | Cross-system SE<->EM connection attempt with QBCP tori; later reached a closed arc #538-#544-#611-#612-#617-#618-#619-#620-#626-#646 whose final negative is a ~166,016 km closure floor; "not a proof of non-existence". | Partly. #538 / #544 / #533 chain used QBCP + BCR4BP seed (exposed, OUTSTANDING.md:912; scripts/run_538_qbcp_cycler.py:90,92,95,99,106). The later methods #611-#646 are listed in the "closed arc" text (OUTSTANDING.md:4729-4750); whether #619/#620/#626/#646 integrate core/qbcp.py could not be established from scripts/ (no script outside run_522/533/538/analyze_593 imports qbcp), so they may have used a different model; could not tell which. | OUTSTANDING.md:4729-4750; OUTSTANDING.md:888,912,1076; scripts/run_538_qbcp_cycler.py:90-106 |
| #544 | Found EM-L2 torus never converges in QBCP (residual 3.4) while SE-L2 converges (3.1e-5); diagnosed BCR4BP mu_sun-continuation seed as degraded (residual 0.139 -> 0.337); proposed per-step convergence gate, finer steps, more Fourier modes; tried the alpha_6 fix on 2026-07-10 and reverted it. | Yes. The diagnosis uses both bcr4bp_torus_residual (wrong Sun sense) and the defective QBCP. #892 withdraws its explanation of the POL1 gap ("two different coefficient sets") - the gap was the parity/scaling defects. | OUTSTANDING.md:7424-7436 (original entry), 872, 909, 1076-1077; docs/notes/2026-07-14-593-qbcp-alpha6-impact-scoping.md:68-76 |
| #593 | Scoped whether the #592 alpha_6 fix changes past QBCP conclusions: |alpha_6-1| up to 0.84%; SE-L2 torus mean state moves by 0.194 between buggy and "fixed"; independent multi-shooting L1-substitute: distance to POL1 4.14e-2 -> 1.81e-2; closed with "#592 stays applied". | Yes, and withdrawn in part. The "fixed" model was itself still defective (#892: wrong parities; #592 scaling "was ours"); OUTSTANDING.md:872 states the #593 explanation of the POL1 gap is withdrawn. The scaling magnitude (0.84%) is arithmetic on the series and independent. The L1 substitute used the bicircular model only as a starting guess (OUTSTANDING.md:968-969, "end result not exposed" to #891). | docs/notes/2026-07-14-593-qbcp-alpha6-impact-scoping.md:8-22,63-90; OUTSTANDING.md:872,888,909-912,968-969; scripts/analyze_593_qbcp_l1_substitute_reconciliation.py:32-33 |

## (iii) scripts/ that use the models (directly or via genome/search modules) AND write results

Writers (model use + file written):

| Script | Model use | Writes | Output exists on disk |
|---|---|---|---|
| scripts/run_303_bcr4bp_l1_continuation.py | :37 core.bcr4bp, :39 genome.bcr4bp_continuation, :42 genome.bcr4bp_genome | :56 OUT_PATH data/bcr4bp_l1_family_303.jsonl; written :270-274 | yes (51 lines) |
| scripts/run_303_catalogue_validation_probe.py | no direct import; reads the 303 family file (:43 FAMILY_PATH), i.e. model output (indirect, downstream) | :45 OUT_PATH data/bcr4bp_validation_bridges_303.jsonl (:228 mkdir) | yes |
| scripts/run_304_bcr4bp_halo_continuation.py | :39 core.bcr4bp, :40 bcr4bp_continuation, :43 bcr4bp_genome | :53 OUT_PATH data/bcr4bp_halo_family_304.jsonl; written :303-307 | yes |
| scripts/run_304_catalogue_validation_probe.py | indirect: reads :52 data/bcr4bp_halo_family_304.jsonl | :54 OUT_PATH data/bcr4bp_halo_validation_bridges_304.jsonl (:264 mkdir) | yes |
| scripts/scan_313_sun_jupiter_moons.py | :75 core.bcr4bp, :77 bcr4bp_continuation, :78 bcr4bp_genome | :453-462 data/scan_313_sun_jupiter_europa.jsonl, scan_313_sun_jupiter_io.jsonl (mars files data/scan_313_mars_*.jsonl also exist; their writer was not traced, could not tell) | yes |
| scripts/scan_334_bcr4bp_system_swap.py | :80 core.bcr4bp, :82 bcr4bp_continuation, :83 bcr4bp_genome, :92 bcr4bp_systems | :103 OUT_PATH data/scan_334_bcr4bp_system_swap.jsonl; written :666-671 | yes |
| scripts/run_538_qbcp_cycler.py | :90 core.bcr4bp, :92 core.qbcp, :95 genome.bcr4bp_torus, :99 genome.qbcp_torus, :106 search.variational_qbcp_torus | :188 data/runlogs/run_538_qbcp_cycler.jsonl, appended :214-216 | yes (runlogs) |
| scripts/screen_884_sun_forced_em_cyclers.py | :39 search.sun_forced_periodic_884 (imports core.bcr4bp at :78 of that module) | :41 OUT = data/found/884_sun_forced_em_cyclers; json written :155-157, :404-437 (controls/, families/). The directory already carries INVALID_WRONG_SUN_SENSE.md | yes |
| scripts/_apply_378_empty_region.py | no model import (only data.empty_regions); it records the #378 verdict whose method used core.bcr4bp | :108 append_empty_region(DEFAULT_EMPTY_REGIONS_PATH, ...) into data/empty_regions.jsonl (JL39) | yes |

Model-using scripts that do NOT write result files (stdout only, or write nothing the survey could find; each grep for open/write/json.dump/save/data paths found nothing):
- scripts/run_305_bcr4bp_gauntlet_probe.py (:29-34 core.bcr4bp, validation v0..v3_bcr4bp; header :13 "REPORTING ONLY ... does NOT write the catalogue"; reads data/bcr4bp_l1_family_303.jsonl, halo_304 :135-136; prints JSON :132,:150)
- scripts/run_522_coherent_connection.py (:17 core.bcr4bp, :21 genome.bcr4bp_torus; prints at :402-405; calls preflight_search at :155 which is a gate, not a result writer - not read in detail)
- scripts/run_533_qbcp_connection.py (:20,22,25,30; prints; preflight_search :218)
- scripts/search_coherent_connections.py (:18,20; no write found)
- scripts/analyze_593_qbcp_l1_substitute_reconciliation.py (:32 core.bcr4bp, :33 core.qbcp; no write found)
- scripts/spike_378_bct_apoapsis.py (:35 core.bcr4bp; prints only). Note its numbers fed JL39 only through the hand-written entry.
- scripts/spike_412_bcr4bp_reach.py (:21-33; prints only, :68-96). Numbers fed JL60 through the entry.

Total scripts importing either model or a module built on them: 16 (matches OUTSTANDING.md:967 "16 scripts"):
analyze_593, _apply_378 (indirect only; no model import), run_303 x2, run_304 x2, run_305, run_522, run_533, run_538, scan_313, scan_334, screen_884, search_coherent_connections, spike_378, spike_412.
Strictly by import: 15 import a model module; `_apply_378` does not but records a model-dependent verdict.

Not exposed (checked, for completeness): screen_694/695/696/701/703/716/882, run_704/705/729, verify_721/724 import variational_ccr4bp_torus / variational_crnbp_torus only (see importer map above).

## Gaps / could not tell
- Which script under scripts/ launched the #378 44-point grid sweep (module search/cislunar_bct_search.py is the producer; no script in scripts/ imports it).
- Writer of data/scan_313_mars_phobos.jsonl / scan_313_mars_deimos.jsonl (scan_313_sun_jupiter_moons.py writes only the Jupiter files at :453-454).
- Whether #619/#620/#626/#646 (end of the #538 arc) integrate core/qbcp.py; no script outside the list above imports it, so they likely use a different model, but the module they call was not identified.
- Whether #334's fitted exponent k and the #412 structural conclusion survive the sense correction (needs a rerun; not done here).
- search/bvp_integral.py:53 imports core.bcr4bp; nothing under scripts/ imports bvp_integral and I did not trace tests/ or other src users.
- search/literature_check.py contains the strings bcr4bp/qbcp; not read, treated as an anchor table rather than a propagator.
