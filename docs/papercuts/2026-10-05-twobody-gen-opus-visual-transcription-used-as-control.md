# A visually transcribed source table was used as control values with no physics self-check

- Date: 2026-10-05
- Agent: twobody-gen-opus
- Seen before: yes. The project rule "no circular goldens" covers where expected values come from. It does not require checking that a transcription is internally consistent.

What happened: `data/sources/hollister-menning-1970-table3.yaml`, read by eye from an image scan, had 30 cells that differ from the print (mostly 3/5 confusions in turn angles). The #942 control was run on it first and failed. Part of that failure traced to read errors that the periapsis formula r_p = mu/V^2 (1/sin(theta/2) - 1) flags at once: the YAML's own V_r and theta disagreed with its own Rmin by up to a factor of 2.
Workaround: a text-layer copy plus a 400-dpi recheck (docs/notes/2026-10-05-hollister-menning-1970-table3-recheck.md); the YAML was corrected in 259efc0d.
Suggested fix: any transcribed source table with redundant columns gets a self-consistency test (a physics identity between its own columns) when it is added, before it is used as expected values.
