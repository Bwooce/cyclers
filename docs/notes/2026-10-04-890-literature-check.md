# #890 literature check: is the Titania-Oberon two-moon periodic flyby orbit already published?

Date 2026-10-04. Scope: the object in `docs/notes/2026-10-04-890-titania-oberon-candidate.md`, a symmetric
periodic orbit of the planar concentric circular restricted four-body problem (Uranus, Titania, Oberon,
both moons at physical mass), period 5 Titania-Oberon synodic periods (123.16 d), one Titania flyby
(1,977 km altitude) and one Oberon flyby (1,364 km altitude) per cycle, reached by continuing a
patched-conic Titania-Oberon-Titania closure in the moon mass. A model object, not a trajectory of the real
system. Assumption throughout: known until I fail to find it.

Evidence grades used below: FULL TEXT (a held paper read in this session, quoted), HELD DIGEST (the
project's own digest of a full text read earlier, not re-read by me), ABSTRACT (abstract or search-result
snippet only), INFERENCE.

## 1. Verdict

Same class as prior work, this body set not found, with three caveats that the coordinator and owner need.

1. No paper found that prints a periodic orbit, or any repeating sequence, that alternates Titania and
   Oberon flybys, in any model. That is a negative over 21 searches (section 4) plus the held corpus, so
   it is conditional on that search and is necessary, not sufficient.
2. The ARCHITECTURE (a free-return cycler between two moons of one planet) is published (Russell and
   Strange 2007/2009) but for Jupiter and Saturn only; the held full text does not treat Uranus.
3. The exact BODY SET and MODEL (Uranus, Titania, Oberon in the planar concentric circular restricted
   four-body problem) WAS computed by Kumar and Anderson (AAS 24-288), with Oberon as the base moon and
   Titania as the perturber. They computed resonant orbits and tori, not flyby cyclers. Under the owner's
   2026-10-03 ruling on "system" in spec 16.4 (ii), this matters (section 5).
4. The THEORY class is published and old: periodic orbits that shadow chains of collision arcs ("second
   species", Poincare; Bolotin and MacKay; Marco and Niederman; Font, Nunes and Simo), including a 3-centre
   version. Their results are for the singular small-mass limit, which is exactly the regime in which the
   #890 note says the patched conic lives (section 2.9 of that note). I have only abstracts and search
   snippets for these, not the full texts.

## 2. What was verified from a full text (held PDFs, read here)

All held PDFs are in the private paper corpus; filenames below are the corpus filenames only.

**Canales, Howell and Fantino 2021** (CMDA 133:36, DOI 10.1007/s10569-021-10031-x; filed in the private
paper corpus as canales-howell-fantino-2021-...-cmda-133-36-...-published.pdf; text extracted and grepped).
Abstract (also fetched from arXiv 2110.03683): "connections between the periodic orbits of such two
different moons are achieved ... Case studies are presented for the Jovian and Uranian systems." Section
4.4 (full text): "Consider a s/c located in an L2 northern halo orbit in the Uranus-Titania (U-T) system
... [arriving at an] L1 southern halo orbit (J Ca = 3.003) in the Uranus-Oberon (U-O) system". The
Fig. 32 caption reads "Transfer from an L2 northern halo orbit of the Uranus-Titania system, to an L1
southern halo orbit of the Uranus-Oberon system ...: Delta-v tot = 45.7 m/s and t tot = 28.44 days".
Contains: a one-way, single-impulse-class transfer between halo orbits of Titania and of Oberon, built
from two coupled three-body problems and a patched analytic connection, then moved to an ephemeris.
Does not contain: any flyby of either moon as a design element, any periodic or repeating Titania-Oberon
sequence, any four-body periodic orbit. The #890 note's description ("a one-way Titania to Oberon
transfer") is accurate.

**Russell and Strange 2009** (JGCD 32(1):143-157, DOI 10.2514/1.36610; filed in the private paper corpus as
russell-strange-2009-cycler-trajectories-planetary-moon-systems-JGCD-32-...pdf). Abstract (verbatim,
fragment): "Previously applied enumerative cycler search and optimization techniques are generalized and
specifically implemented in the Jovian and Saturnian moon systems. Overall, hundreds of ideal model ... are
found". Table 1, "Planetary moon ideal models considered": Jupiter: Ganymede to Io, Ganymede to Europa,
Ganymede to Callisto, Europa to Ganymede; Saturn: Titan to Enceladus. Text: "Titan is capable of providing
the gravity-assists as it is the largest moon in the Saturn system by almost 2 orders of magnitude.
Furthermore, the small mass of Enceladus validates the massless assumption of the target body in the ideal
model". Contains: free-return cyclers with ONE working flyby body and a massless target, found in a
circular coplanar patched-conic model, then transitioned to a fully integrated ephemeris model (the
abstract says near-ballistic cyclers result). Does not contain: Uranus (a grep of the extracted text finds
no Uranian system treated), any restricted-four-body periodic orbit, any mass continuation. The #890
object has BOTH moons working (both at physical mass, both turning the trajectory by 66 and 70 degrees),
which is outside this architecture's one-working-node assumption.

**Strange, Russell and Buffington 2007** (AAS 07-277, held full text). Line 41 of the extracted text:
the Tisserand-criterion method is "extremely useful for the heliocentric case [6] as well as for Jovian [7]
and Uranian [8] tour design", where [8] is Heaton and Longuski 2003. Resonance hopping between two
circular coplanar moons; no Uranian cycler.

**Strange, Landau and Longuski 2013** (AAS 13-801, held full text; Uranian inclination-reduction
sequence). "We see that a Titania 3:1 resonance (26 days) is very close to an Oberon 2:1 resonance (27
days)." This is a Titania-Oberon Tisserand-graph near-coincidence at a roughly 26 d orbit period, used for
one-way inclination reduction, not a cycler. The #890 legs (about 11 d orbit period, 5.5 revolutions, 61.58
d) are a different neighbourhood.

**Heaton and Longuski 2003** (JSR 40(4):591-596, DOI 10.2514/2.3981; held, OCR text). A one-shot 811-day
three-phase Galileo-style tour with more than 40 flybys of Ariel, Umbriel, Titania and Oberon; Oberon
flybys listed in a flyby table (OCR is poor). It is a tour that ends, not a periodic orbit.

## 3. Held digests (read in this session, full text not re-read)

**Kumar and Anderson, AAS 24-288, 2024** (private corpus file anderson-kumar-2024-oberon-mmr-unstable-orbit-
survey-aas-24-288.pdf; digest `docs/notes/2026-07-27-728-anderson-kumar-2024-oberon-mmr-survey-digest.md`,
HELD DIGEST of a full text with page numbers). Uranus-Oberon PCRTBP resonant orbits and heteroclinics,
then "Uranus-Oberon-Titania concentric circular restricted 4-body problem (CCR4BP)" with Titania at
mu3 = 3.91677e-5, continuing Oberon resonant orbits to tori or secondary-resonant periodic orbits.
Conclusion (p.19, as quoted in the digest): "In this study, we do not yet look at the heteroclinics between
mean motion resonances in the CCR4BP", and "there are three other large moons of Uranus - Titania, Umbriel,
and Ariel - for which this study should be repeated". Contains: the same four-body model, the same body
set, periodic orbits that are Oberon-resonant (a spacecraft orbiting Uranus, not flying by two moons in
alternation). Does not contain: Titania flybys as design elements, a Titania-Oberon flyby cycle. The
digest also records that the 4:3 Oberon family intersects Titania's orbit, so no torus continuation exists
there; the paper does not compute those orbits.

**Lynam and Longuski 2011** (Acta Astronautica 69; digest `2026-06-30-digest-lynam-longuski-2011-laplace-
resonant-triple-cyclers.md`) and **Hernandez, Jones and Jesick 2017** (the held file is AAS 17-608, the
Io-Europa-Ganymede one-class triple cyclers; the task brief says AAS 17-462, I did not resolve which number
is correct; the held digest is `2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md`): Jovian
Laplace-resonant triple cyclers in patched-conic models. A grep of the extracted Lynam-Longuski text finds
only generic mentions of Uranus-mass exoplanets and rings, no Uranian design. No periodic-orbit
refinement in a restricted four-body model.

**Anderson and Lo 2010/2011** (digest `2026-07-28-745-anderson-lo-2010-2011-resonant-flyby-digest.md`):
a ballistic Jupiter-Europa resonance-cycling flyby trajectory shown to shadow a homoclinic/heteroclinic
connection in the CR3BP. One moon, no four-body periodic orbit.

**Kumar, Anderson, de la Llave et al. (Jovian CCR4BP series)** (digests of 2021, 2023): Jupiter-Europa and
Jupiter-Ganymede resonant tori and transfers in the CCR4BP. Same model family as #890; Jovian; orbits are
spacecraft-resonant, not alternating-flyby cycles.

## 4. Live search, every query

WebSearch unless marked. "Empty" means nothing on point.

1. "Titania Oberon cycler trajectory Uranus moons periodic flyby": Canales et al. 2110.03683, AAS 24-288,
   dynamics papers. No cycler. Empty for the object.
2. "Russell Strange Cycler trajectories in planetary moon systems Uranus": the 2009 paper; no Uranus
   mention in the snippets (confirmed absent in the held full text).
3. "Uranus orbiter tour Titania Oberon resonance hopping gravity assist V-infinity leveraging": UOP trajectory
   options, McAdams 2011, AAS 24-288. Tour design, one-way; "Titania and/or Oberon" used for inclination
   cranking. Empty for periodic.
4. "periodic orbits restricted four-body problem two moons flyby sequence continuation patched conic
   cycler": CCR4BP Jovian (Kumar), Hill restricted 4-body (arXiv 2402.19181), a "Periodic orbits in the
   restricted four-body problem" paper (ScienceDirect, 1986, Lagrange-configuration), Koon-Lo-Marsden-Ross
   book. None on flyby-chain periodic orbits.
5. "Bolotin MacKay periodic orbits near collision chains second species restricted three-body problem":
   Bolotin and MacKay, "Nonplanar second species periodic and chaotic trajectories for the circular
   restricted three-body problem" (CMDA 94, 2006, DOI 10.1007/s10569-006-9006-0); Bolotin, elliptic second
   species (CMDA, DOI 10.1007/s10569-005-2172-7); Barrabes-Cors-Pinyol-Soler, Guardia et al. on oscillatory
   and collision orbits. ABSTRACT level.
6. "Henon second species ... Font Nunes Simo": Font, Nunes, Simo, "Consecutive quasi-collisions in the planar
   circular RTBP" (Nonlinearity 15(1):115-142, 2002) and "A numerical study of the orbits of second species of
   the planar circular RTBP" (CMDA 103(2):143-162, 2009), mass ratio treated as perturbation up to 1e-3;
   "A general timing condition for consecutive collision orbits ... elliptic restricted problem" (Springer).
   ABSTRACT level.
7. "Lynam Longuski multiple-satellite-aided capture ...": Jovian double, triple, quadruple capture;
   CGI broad-search papers. Jovian only.
8. "Hernandez Jones Jesick AAS 17-462 Galilean triple cyclers": Io-Europa-Ganymede triple cyclers, 2017;
   Liang et al. CGE triple cyclers. Jovian only.
9. "Strange Russell Buffington Mapping the V-infinity globe": confirms Saturn/Jupiter multi-moon resonance
   hopping; no Uranus.
10. "Uranian satellites ballistic cycler between Titania and Oberon free-return": returned this project's own
    public repository (Bwooce/cyclers) and AAS 24-288. See section 7 on the repository hit.
11. "periodic orbit flybys two moons same planet concentric circular restricted four-body problem Galilean
    moons": Kumar et al. AAS 21-651 (arXiv 2109.14815), 2109.14800, 2309.06073. Resonant tori, not
    alternating-flyby cycles.
12. "Uranus flagship tour design 2025 Titania Oberon flybys repeated moon tour trajectory paper": Uranus
    Orbiter and Probe concept update (PSJ ae680c), McAdams 2011, Landau/Strange. Tours; "repeated flybys of
    Titania" for inclination reduction; no cycler.
13. "Titania Oberon cycler Uranus moons trajectory arXiv 2025 2026": moon discovery news, MUSE concept. Empty.
14. "Heaton Longuski Uranian satellites tour Oberon Titania resonance Galileo-style repeated flybys":
    confirms the one-shot tour; Titania and Oberon as the Ganymede-Callisto analogues.
15. "Marco Niederman periodic orbits consecutive close approaches ... shadowing collision orbits": Marco and
    Niederman, "Sur la construction des solutions de seconde espece dans le probleme plan restreint des trois
    corps" (Ann. IHP Phys. Theor. 62:211-249, 1995); Bolotin and MacKay symbolic dynamics; Dimare 3-centre
    (arXiv 0911.3557). ABSTRACT level.
16. WebFetch arxiv.org/abs/0911.3557: Dimare, "Chaotic quasi-collision trajectories in the 3-centre problem":
    "uniformly hyperbolic invariant sets of periodic and chaotic almost collision orbits" for a small third
    centre, via Bolotin-MacKay theory. Abstract only. Relevant to a three-centre (two fixed primaries plus
    a third) setting, not a rotating two-moon problem.
17. "Bolotin Negrini regularization and topological entropy spatial n-center problem": Bolotin and Negrini,
    Ergodic Theory Dynam. Systems 21(2):383-399, 2001; a snippet states that in the plane 3-body problem with
    2 small masses, second species solutions "shadow chains of collision orbits of 2 uncoupled Kepler
    problems". I did not see the source of that sentence; treat as a lead, SNIPPET ONLY.
18. "concentric circular restricted four-body Uranus Titania Oberon periodic orbit flyby": only 2110.03683 and
    AAS 24-288. Empty for flyby periodic orbits.
19. "Hernandez Jones Jesick Jones Longuski Uranus moon cycler ... triple cycler Uranian": nothing Uranian.
20. "Europa Ganymede cycler converted to periodic orbit concentric circular restricted four-body problem": Kumar
    et al. series; "Planetary moon cycler trajectories" (Russell and Strange, AAS/AIAA Space Flight Mechanics
    Meeting, Sedona, Feb 2007; JPL TRS 2014/40318, the conference precursor of the 2009 paper); a Journal of
    the Astronautical Sciences paper titled "A Continuation Method for Converting Trajectories from Patched
    Conics to Full Gravity Models" (DOI 10.1007/s40295-014-0017-x; I saw only the title, the publisher page
    redirected to a login, authors and abstract not seen). That title is the closest method-level hit for the
    patched-conic-to-full-model continuation step; whether it continues in a mass parameter or in a
    perturbation ladder is unknown to me.
21. "Bwooce cyclers Titania-Oberon-Titania quasi-cycler Uranian symmetric closure": this project's own
    repositories only (section 7).
22. WebFetch arxiv.org/abs/2509.03655 (Kumar, multi-shooting parameterization, J. Nonlinear Sci. 2026,
    DOI 10.1007/s00332-026-10276-6): abstract has no Titania, Oberon, flyby or cycler content.
23. WebFetch arxiv.org/abs/2110.03683: abstract as quoted in section 2.
24. WebFetch researchgate Russell-Strange 2007 page: HTTP 403, not read. WebFetch of the Springer JAS page:
    redirect, not read.

That is 21 distinct WebSearch queries (list items 20 covers two) plus fetch attempts: three arXiv abstract
fetches that worked (2110.03683, 2509.03655, 0911.3557), one irrelevant fetch of arXiv 2111.11858 that I
discarded, and two blocked (ResearchGate 403, Springer redirect).

## 5. Which spec 16.4 case I think applies (the label is the coordinator's and owner's call)

Literal collision first: no published orbit, family, or tabulated member found that this orbit belongs to,
so the literal-collision test does not fire on what I found. Then:

- Case (ii), known architecture at a never-treated system, is the nearest fit but is NOT clean. Russell
  and Strange's computed sets (Table 1: Ganymede to Io/Europa/Callisto, Europa to Ganymede, Titan to
  Enceladus) do not include Uranus, so the architecture is "applied at a body set the source never
  computed". The ruling reads "system" as "the BODY SET the source actually computed" and the source here is
  R-S. But the #890 object is not in R-S's architecture (two working nodes, physical masses, a four-body
  periodic orbit rather than a patched-conic free return), so "applying Russell and Strange's architecture
  unchanged" is not the claim wording that is true. The claim that is true is closer to "a periodic orbit of
  the Uranus-Titania-Oberon concentric circular restricted four-body problem that continues a patched-conic
  closure".
- The body set and model were computed by Kumar and Anderson (AAS 24-288). If the owner's test is
  "did a source compute this body set in this model", the answer is yes (for resonant orbits and tori); if
  it is "did a source compute this orbit class here", no. The ruling says a set the source "did not compute"
  is never-treated "unless the source specifically claims to have searched it". Kumar and Anderson did not
  claim to have searched for flyby cycles, and say heteroclinics in the CCR4BP are future work. I read
  that as not blocking (ii), but it is the sentence in the ruling most worth re-reading against this case.
- Case (iii), theorem-generic object, may apply to the existence of the orbit: second-species theory
  (Bolotin and MacKay; Marco and Niederman; Font, Nunes and Simo) gives periodic orbits shadowing collision
  chains for small mass ratio. I could not verify from a full text that a theorem covers two small moons
  of one primary with distinct periods in a rotating frame (the Bolotin-Negrini snippet and Dimare's 3-centre
  result point that way). If the theory does cover it, the policy wording is "first computed", never "first
  predicted", and the parent family is the collision-chain family.
- My own recommendation: treat as case (ii) with the wording above, state in `notes` that the same model and
  body set were treated for Oberon-resonant orbits by Kumar and Anderson, and do not describe the object as
  new in kind, because the singular-limit theory is a published class. Both (ii) and (iii) leave it
  `candidate-novel` under the policy and neither is a novelty-proof; the project's own note already says
  it is a model object and not a catalogue row.

## 6. The project's gate: `check_literature`

Called with `primary="Uranus"`, `sequence=("Titania","Oberon","Titania")`, `period_k=5`,
`vinf_per_encounter_kms=(0.274732, 0.268904, 0.274732)`, `n_rev=(5,5)` and `search=offline_corpus_search`
(the corpus-anchor search; the live-web path needs an injected WebSearch). Script kept in the session
scratchpad, not in the repo. Output:

```
status='published'  confidence=0.95  doi='10.2514/2.3981'
matched_url='https://doi.org/10.2514/2.3981'
query_trail=['Titania-Oberon cycler 5 synodic trajectory']
citation="Heaton & Longuski, 'The Feasibility of a Galileo-Style Tour of the Uranian Satellites,' J. Spacecraft
 & Rockets 40(4):591-596 (2003); AIAA 2001-3859 / NTRS 20020021945. Foundational one-shot Uranus moon-tour
 anchor: 811-day three-phase tour with 40+ flybys across Miranda/Ariel/Umbriel/Titania/Oberon.
 NOT a periodic cycler."
notes='Structural fingerprint matched a published cycler -- treat as a rediscovery; NOT novelty-claimable.'
is_novelty_claimable(...) = False
```

Reading: the gate is not "clear". It returned `published` against the Heaton and Longuski anchor, and the
anchor's own text says "NOT a periodic cycler". I read that as a body-set false positive (Uranus plus a
moon sequence overlap), not as evidence the orbit is known; the hit is the wrong class for the reason in
section 2. Only one query was generated, so the gate did not exercise the live web. The gate is a filter,
not a proof, in either direction.

## 7. Anything the literature contradicts in the #890 note

Nothing contradicted. Points to correct or add:

- The note says the Canales-Howell-Fantino treatment is "a one-way Titania to Oberon transfer": confirmed
  (halo to halo, 28.44 d, 45.7 m/s). It may add that the paper's moon models are coupled three-body
  problems, not a four-body model.
- The note does not cite Kumar and Anderson AAS 24-288, the one paper that computes periodic orbits in this
  exact four-body model with this body set. It should, with the base/perturber role difference (Oberon base,
  Titania perturber there; Titania base, Oberon perturber in #890). The note's remark that the 4:3
  Oberon-family crosses Titania's orbit is not in the note and is not needed.
- The note's "singular limit" account of the patched conic (section 1.3, 2.9) and the observed sub-linear
  growth of periapsis distance with mass match the collision-chain picture from second-species theory
  qualitatively (arcs approaching collision as the mass goes to zero). That is an inference from abstracts,
  and the theory would be the natural independent check on the continuation behaviour, not a contradiction.
- Search-engine summaries returned this project's own public repository (Bwooce/cyclers, and
  Bwooce/cyclers.space) describing a Titania-Oberon quasi-cycler as "discovered and validated". That is
  project output, not prior art, and I did not fetch the pages. Given the 2026-10-04 withdrawal of the six
  Uranian quasi-cycler rows (`#888`), someone should check that the public README text is not stale.
- The task brief gives the Hernandez, Jones and Jesick paper as AAS 17-462; the held file and digest are
  AAS 17-608. Unresolved; both are Jovian.

## 8. Not held by the project; open-access or findable

Held and sufficient: Canales-Howell-Fantino 2021, Russell-Strange 2009, Strange-Russell-Buffington 2007,
Strange-Landau-Longuski 2013, Heaton-Longuski 2003, AAS 24-288, Lynam-Longuski 2011, Hernandez 2017.

Not held (worth filing, in priority order):

1. Bolotin and MacKay, "Nonplanar second species periodic and chaotic trajectories for the circular
   restricted three-body problem", CMDA 94 (2006), DOI 10.1007/s10569-006-9006-0. Decides whether the theory
   case (iii) applies and what "second species" gives for two moons.
2. Font, Nunes and Simo, "Consecutive quasi-collisions in the planar circular RTBP", Nonlinearity 15(1):
   115-142 (2002), and "A numerical study of the orbits of second species of the planar circular RTBP",
   CMDA 103(2):143-162 (2009). Numerical computation of the orbit class with mass ratio as the
   continuation parameter up to 1e-3: the nearest published analogue of the #890 continuation, for one
   moon.
3. Bolotin and Negrini, "Regularization and topological entropy for the spatial n-center problem",
   Ergodic Theory Dynam. Systems 21(2):383-399 (2001), and Dimare, arXiv 0911.3557 (open access): the
   n-centre and 3-centre quasi-collision results.
4. Marco and Niederman, Ann. IHP Phys. Theor. 62:211-249 (1995), in French.
5. Russell and Strange 2007, "Planetary moon cycler trajectories", AAS/AIAA Space Flight Mechanics Meeting
   (Sedona); copy at JPL Technical Reports Server, handle https://trs.jpl.nasa.gov/handle/2014/40318. The
   conference precursor of the 2009 paper; it should be checked for a Uranus case.
6. Kumar, arXiv 2509.03655, DOI 10.1007/s00332-026-10276-6 (open access on arXiv). Abstract has no Titania
   content; it is the method companion to AAS 24-288.
7. The JAS "A Continuation Method for Converting Trajectories from Patched Conics to Full Gravity Models",
   DOI 10.1007/s40295-014-0017-x. Title only seen; read before claiming the patched-conic-to-periodic-orbit
   continuation is not a published method.

## 9. Limits of this check

- 21 searches plus the held corpus is a search of indexed text; conference proceedings and theses are
  under-indexed and the second-species literature is mostly mathematics, not astrodynamics, so a
  mission-design paper that mentions a Titania-Oberon cycle in passing would not necessarily be found.
- The second-species and n-centre papers were seen only as titles, abstracts and search-result snippets.
  Whether a published theorem or numerical study already covers two moons of one planet is the main
  unresolved question and needs the papers in section 8, items 1 to 3.
- Several held digests (Kumar and Anderson, Lynam-Longuski, Anderson-Lo, Kumar Jovian series) were read as
  digests, not re-read from the PDFs.
