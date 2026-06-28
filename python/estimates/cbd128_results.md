# CBD Parameter Search Results (128-bit security)

Cost model: BDGL16

## n=636, Xs=Binary

| q     | min eta | security (bits) |
|-------|--------:|----------------:|
| 2^12  |       4 |           131.9 |
| 2^14  |       4 |           131.7 |
| 2^15  |       8 |           128.0 |
| 2^16  |      31 |           128.0 |
| 2^18  |     N/A |           117.8 |
| 2^20  |     N/A |           106.1 |

## n=1024, Xs=UniformMod(3)

| q     | min eta | security (bits) |
|-------|--------:|----------------:|
| 2^16  |       4 |           144.7 |
| 2^20  |       4 |           144.7 |
| 2^24  |       4 |           138.2 |
| 2^25  |       4 |           132.6 |
| 2^26  |     N/A |           127.3 |
| 2^27  |     N/A |           122.5 |
| 2^28  |     N/A |           127.3 |
| 2^32  |     N/A |           110.3 |

Notes:
- "min eta" is the smallest CenteredBinomial(eta) parameter achieving >= 128-bit security at that q.
- "N/A" means even eta=64 does not reach 128-bit security.
- For n=636, the dominant attack shifts from arora-gb (small eta) to dual_hybrid (large q).
- For n=1024, eta=4 suffices up to q=2^25. Beyond that, lattice attacks (dual_hybrid) break 128-bit.
