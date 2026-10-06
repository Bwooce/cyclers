# #942 / #943 ge-* and em-* candidates: web prior-art search (2026-10-07, lower priority)

Requested by the lead as optional, because these candidates are not proposed as finds. ge-1..3 are
real-ephemeris negatives (conditional), and em-1..3 do not pass rung (d). See
`2026-10-06-942-943-owner-decision-summary.md` sec. 3. Search by twobody-gen2-opus.

## Method

- OpenAlex title/abstract searches, 11 queries:
  - ge: "ganymede europa cycler", "europa ganymede periodic orbit tour", "ganymede europa resonant
    repeating trajectory", "jovian moon cycler", "double cycler jupiter", "europa ganymede gravity
    assist periodic";
  - em: "earth mars cycler mars gravity assist", "mars swingby periodic orbit earth", "earth mars cycler
    mars flyby", "two planet cycler mars venus earth both flybys", "mars flyby cycler ballistic".
- The held Jovian collision checks done by corpus-file-opus (results note 6.41, 6.44):
  - 21F31 (ISSFD 2024), Cangahuala 2025, Lam 2015, Lam 2018;
  - Buffington 2012, Buffington 2014;
  - Golubev 2014, Grushevskii 2017;
  - Minovitch 1972, Boutonnet 2024, Liang 2024.
- The Earth-Mars forward-citation sweeps of 2026-06-11 and 2026-06-13.

## ge-1, ge-2, ge-3 (Ganymede-Europa, both bend)

- No two-working-body Ganymede-Europa cycler was found in the literature.
- The hits are:
  - Russell & Strange 2007/2009: the one-body GanEur/EurGan rows; the ge cell recovers GanEur#43 and
    EurGan#131 in-run;
  - Liang et al. 2024: CGCEC triple cyclers, three moons;
  - Lynam 2012 (PhD) and Lynam & Longuski 2011: Laplace-resonant Io-Europa-Ganymede triples;
  - Haapala 2018 "Patched periodic orbits" (UT Austin thesis): CR3BP periodic orbits for moon tours,
    another model;
  - Kumar, Anderson & de la Llave 2023: G-E resonant tori, CR4BP, context.
- The Clipper tours have no repeating G-E segment: one pump-down handoff, and the 11-F5 switch-flip
  used once.
- **Verdict: no collision.**

## em-1..em-5 (Earth-Mars, both bend)

- Hits:
  - the Russell-Ocampo and McConaghy-Russell-Longuski families: Earth hosts, Mars massless;
  - Rall 1969 / Rall & Hollister 1971: Mars swing-bys of 2.3-13.6 deg, k = 4-6, no Mars returns;
    already checked in 6.36;
  - Pisarevsky 2008: method; Table 4 class III; Fig. 14 not digitised;
  - Jones 2017 VEM triples;
  - Rogers 2014, McConaghy 2004 and Hughes 2016 PhDs;
  - "An Alternative Humans to Mars Approach: Multiple Mars Flyby Trajectories" (AIAA 2015-4412):
    one-shot.
  None is a two-working-body Earth-Mars cycler with a generic Mars return.
- **One unchecked possible relative, NOT HELD:** Fornari & Pontani 2020, "A global exploration method
  to identify families of cycling Earth-Mars trajectories", Aerotecnica Missili & Spazio, doi
  10.1007/s42496-020-00050-6. It is a global search for cycling E-M families, and its OpenAlex record
  has no abstract. Whether it allows Mars to bend is unknown. It has been gated (Springer) since the
  2026-06-11 forward-citation sweep.
  - Also not held: Pelle et al. 2019, "Earth-Mars cyclers for a sustainable human exploration of Mars",
    Acta Astronautica 154:286 (an architecture trade, low risk).
- **Verdict: no collision found among held and open sources.** Fornari & Pontani 2020 is the one
  source to acquire before any em row could be proposed. No em row is proposed now.

## Acquire (priority order)

1. Fornari, E. & Pontani, M. (2020), Aerotecnica Missili & Spazio, doi 10.1007/s42496-020-00050-6 (em).
2. Pelle, S. et al. (2019), Acta Astronautica 154:286-294 (em, low risk).
