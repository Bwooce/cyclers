# Triage: Rickman, NASA NESC slides — Intro to Orbital Mechanics & Spacecraft Attitudes for Thermal Engineers

**Date:** 2026-06-10
**Source:** Rickman, NASA NESC slides (203 slides, 5 parts)
**Verdict:** NOT USEFUL — do not re-mine

## What it is

Introductory training deck for thermal engineers covering: vector/matrix math, two-body
problem derivation, Kepler's laws, perturbed orbits (J2, node precession, beta angle,
eclipse fractions), transfer orbits / gravity assists / restricted three-body / halo orbits
(brief), and spacecraft attitude transformations (LVLH, Celestial Inertial, Euler sequences).

Primary source throughout is Bate, Mueller & White, *Fundamentals of Astrodynamics*
(Dover, 1971), supplemented by Wikipedia and MathWorld. No independent reference list
slide exists; all citations are inline.

## Check 1 — Golden-test numeric values

None usable. All worked examples (solar flux vs. Kepler's 1st law, geostationary orbit
period, ISS eclipse fraction, beta-angle time-series) are thermal-environment problems
with approximate constants. Slide 3 disclaimer explicitly states the values "should not
be used for design, analysis or mission planning purposes." No trajectory propagation
examples with published input+output pairs appear anywhere in the deck.

## Check 2 — Reference list pointing at relevant sources

No references to cyclers, Earth–Mars trajectories, interplanetary mission design, or
Lambert solvers. The single non-trivial reference is Bate/Mueller/White, which we already
hold. No new pointers.

## Summary

Pure introductory pedagogy aimed at thermal engineers. Nothing overlaps with cycler
trajectory finding, Tisserand graphs, or mission-design numerics. Safe to ignore permanently.
