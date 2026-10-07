# The two Bruno-Varin family-h source tables print the period in different conventions

- Date: 2026-10-07
- Agent: corpus-review2-fable (found by its reader C1 during #971)
- Seen before: no (the digests of the two papers each note their own convention; nobody had compared them)

What happened: `data/sources/bruno-varin-2009-family-h-small-mu-tables.yaml` (SSR 43(1):2, mu = 0, 0.00095,
0.1, 0.2) prints T~ = (1 - mu) T / 2pi although the paper's text says T / 2pi; `data/sources/bruno-varin-2009-family-h-big-mu-tables.yaml`
(SSR 43(2):158, mu = 0.3, 0.4, 0.5) and `varin-bruno-2009-closed-families-tables.yaml` print T~ = T / 2pi.
A DOP853 integration from the printed starts confirmed both: mu = 0.1 Table 4 row 1 reaches its first
x2 = 0 crossing at t/pi = 0.4027 = printed 0.36243 / (1 - mu), while mu = 0.3 Table 1 row 1 reaches it at
t/pi = 0.38873 = printed T~ with no factor. A loader or test that reads the two family-h yamls with one
convention is 25-30 percent wrong in the period at mu = 0.2-0.3 and will reject every row.
Workaround: the small-mu yaml header states its convention; the reader applied it per file by hand.
Suggested fix: add a per-table `period_convention` field (`T_over_2pi` or `one_minus_mu_T_over_2pi`) to
all five Bruno-Varin yamls and make any future loader read it; the same files also print C = -2H + mu
(barycentric Jacobi + mu(1 - mu)), stated in the headers, which the field could carry too.
