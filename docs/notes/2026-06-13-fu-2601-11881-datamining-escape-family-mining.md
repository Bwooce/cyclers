# Fu, Wu, Gong & Shi 2026 — "Data Mining-Based Cislunar Escape-Family Analysis in The Multi-Body Models" (discovery-capability digest)

**Date:** 2026-06-13
**Task:** #250 discovery-capability sweep (head-start candidate 3, Track B).
**Source (free):** S. Fu, D. Wu, S. Gong, P. Shi, School of Astronautics,
Beihang University (BUAA). arXiv:2601.11881 (submitted 17 Jan 2026; v-on-arXiv
title "Data Mining-Based Cislunar Escape-Family Analysis in The Multi-Body
Models"; HTML at arxiv.org/html/2601.11881).
**Access:** abstract + arXiv HTML metadata via WebSearch; `WebFetch` denied
(no figure/table transcription).
**Writeback: NONE** (digest only).

---

## Verdict: BACKGROUND, prioritizer-method-adjacent (Track B) — worth a closer read if/when the discovery daemon needs a clustering layer

> **CONFIRMED from full HTML 2026-06-13 (#251, WebFetch re-enabled).** Full text
> confirms: clustering algorithm is **DBSCAN**; output is one-way *escape*
> trajectories (criterion E>0 + monotonically increasing distance), NOT cyclers;
> Tables 4–6 list 19/21/24 escape families by ΔV/TOF (escapes, not
> catalogue-relevant); cross-model CR3BP→BCR4BP persistence studied (families
> emerge/disappear under solar gravity — concrete examples: Families XX/XXI/XXII–
> XXV emerge, XI/XVI/XVIII/XIX disappear). No cycler-discovery capability; DBSCAN
> family-clustering + persistence is the only transferable idea. Verdict
> unchanged: BACKGROUND (Track-B-conditional).

## What it does

Systematic analysis of **escape trajectories** from a 167 km circular Earth
orbit, in:
1. the Earth-Moon **planar CR3BP**, and
2. the Sun-Earth/Moon **planar bicircular four-body problem** (BCR4BP).

It proposes a method to **identify "escape families"** by combining dynamical
analysis with **data-mining (clustering) techniques** on the escape solution
space. Once families are extracted, it characterizes them (generalized energy,
transfer characteristics) and studies how **solar gravity perturbation** changes
the *count* of escape trajectories and drives the **emergence and disappearance
of escape families** between the CR3BP and BCR4BP models.

## Discovery vs anchor

**Neither directly — the OUTPUT is escape trajectories, not cyclers.** Escape
families are the opposite topology to a repeating cycler (they leave the system;
cyclers return indefinitely). No cycler ICs, no catalogue rows.

The **METHOD is the interesting part for Track B (the prioritizer).** The
discovery spec's Track B asks for a way to point search at non-empty regions and
to cluster a solution space into families. This paper is a concrete instance of:
*sample the solution space → cluster into dynamically-coherent families →
characterize each family by energy/transfer priors*. That is precisely the
"data-mine the solution space for families + initial-state/energy priors"
capability the brief flagged. Two ideas are transferable to our discovery
daemon's dedup/prioritization layer:

1. **Family extraction by clustering**, not by hand-labelled topology — useful
   for auto-grouping daemon candidates before feeding the gauntlet, and for
   defining the negative-registry's "region" boundaries empirically rather than
   by genome parameter box.
2. **Cross-model family persistence (CR3BP → BCR4BP under solar gravity).** The
   "emergence/disappearance of families when SGP is switched on" analysis is the
   same cross-fidelity question our cyclers face (a PCR3BP cycler may or may not
   persist under solar perturbation — cf. the Ross (3,2) family's Leiva-Briozzo
   quasi-bicircular persistence). A clustering-based persistence diagnostic could
   become an automated "does this candidate survive to BCR4BP" screen.

## Reproducible data

Methodology + characterization figures; the abstract/HTML metadata expose no
ICs or family tables we can transcribe (and WebFetch is denied here). The
clustering pipeline (feature definition on escape solutions + the
emergence/disappearance metric) is the reproducible asset, not data rows.

## Method capability

No new cycler genome or topology. As a **prioritizer/clustering method** it is a
candidate *technique* for Track B and for the discovery daemon's
family-dedup/region-definition layer — but it is applied to escapes, so adoption
would mean lifting the clustering approach and re-pointing it at cycler-space,
not using the paper's results.

## Proposed follow-on

**Background, with one conditional action:** if/when the Track-C discovery daemon
needs an automated family-clustering + cross-fidelity-persistence layer (to
group candidates and define empirical empty-region boundaries), pull the full
PDF (user fetch, free) and evaluate this clustering pipeline as a reference
implementation. Until the daemon reaches that stage, do not invest — our current
dedup is catalogue+negative-registry matching, which suffices for the present
seed volumes. Same-group context: Liang et al. (NUAA, not BUAA) is a different
lab; Gong's group (BUAA) is worth a forward-citation watch for cislunar
data-mining follow-ons.
