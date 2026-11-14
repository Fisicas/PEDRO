# VEDA 4.99e on our 614-atom ONIOM model (13 November 2026)

`veda_session_skra_v01.log` is VEDA's own session log (`skra.log`, local folder names removed), in order:
`block14.fmu` (raw 14-atom benzaldehyde block; *Create DD2*, *Calculate*), `near480.fmu` (480-atom block: read, then
`CAUTION: shortage of coordinates 451 instead of 1434`, `error in DD2, coordinate No. 452`), `near481.fmu` and
`full614.fmu` (only `Input file:` is written; VEDA stops without a message). VEDA was restarted before each of the
last two inputs.

`limits_results_v01.json`: PEDRO vs VEDA on the 14-atom block and on UniMoVib's rev-GSVA benzaldehyde Hessian with the
same coordinate set, and PEDRO's C=O stretch TED in cases F and P of the full model. Inputs were written with
`pedro.io.veda.write_fmu` from `validation/own/oniom_benzaldehyde_water/` (see docs/WHY_PEDRO.md).
