# Software and benchmark archives

PEDRO's software and owned ONIOM benchmark are separate research objects. Archive the software through the
Zenodo–GitHub integration, and upload the benchmark as a separate dataset. Link their records after their
identifiers are available. Collaborator-owned research logs are excluded from both archives.

## Software release

`CITATION.cff` is the source of citation metadata for GitHub and Zenodo: Michael Raymond Nelson, MIT license,
software version, title, and description. A separate `.zenodo.json` is unnecessary for the current metadata;
if one is introduced later, Zenodo gives it precedence over `CITATION.cff`.

1. Enable `Fisicas/PEDRO` in the [Zenodo GitHub settings](https://zenodo.org/account/settings/github/).
2. Set `date-released` in `CITATION.cff` and date the changelog on the actual publication day. Commit and push
   the reviewed metadata and code, and require the resulting CI run to pass.
3. Create the GitHub release `v0.2.0b1` from that reviewed commit, titled **PEDRO 0.2.0b1 — public beta**.
   Identify it as a pre-release. Use the [release notes](releases/0.2.0b1.md).
4. Check the release in Zenodo's repository integration page. Verify the archived version, creator, MIT license,
   description, and source archive before announcing the DOI. Check the integration error details if it fails.
5. Add the actual version DOI to `CITATION.cff` and the README citation. Use Zenodo's all-versions DOI badge for
   the project README, and label the two identifiers clearly. Commit these follow-up metadata changes without
   moving the published tag. The initial archive can precede the commit that records its assigned DOI.

For a scientific analysis, cite the **specific software version DOI**, and record the version and analysis
options. The all-versions DOI identifies the software project across releases. A moving `main` branch is not
a substitute for a versioned archive. For unreleased development code, record its full Git commit as well.

There is no software DOI in this document until a record has actually been published. Do not create a second
manual software deposit for the same release while the GitHub integration is processing it.

## Owned ONIOM benchmark dataset

Suggested record title: **PEDRO ONIOM benchmark: benzaldehyde in 200 UFF waters (614 atoms)**.
Creator: **Michael Raymond Nelson**, [ORCID 0009-0004-9809-9842](https://orcid.org/0009-0004-9809-9842).
Resource type: **Dataset**. Dataset version: **1.0** (independent of the software version).
License: **CC BY 4.0** for the owned data, while included software scripts retain MIT licensing. The dataset
contains only the model generated for PEDRO and its derived analysis; collaborator and Q-drive research data
are excluded.

Include the original Gaussian input and complete output log, deterministic model generator, SHA-256 checksums,
PEDRO analysis options and validation record, and the PED/TED tables and unrounded arrays. Preserve original
files and identify derived files explicitly. The raw log must have SHA-256
`40e765cdd93e90b1ffe7dcb3f1a25f5f204fadd7cad529fac04bdb576c3fd131` (42,513,569 bytes).

The model is mechanically embedded `ONIOM(B3LYP/6-31G(d):UFF)`, with benzaldehyde atoms 1–14 and 200 waters.
It is a numerical scale and extraction benchmark, not a physical solvation model or a validated protein/cofactor
assignment. Document case F, unscaled frequencies, signed TED, diagonal versus row-partition PED, and the
approximate environment identity. The [benchmark description](../validation/own/oniom_benzaldehyde_water/README.md)
contains the measured results and limitations.

After publication, add the actual dataset DOI to the benchmark README and validation documentation. On Zenodo,
relate the dataset to the software as **Is supplement to**, and the software to the dataset as **Is supplemented
by**. Cite both records when using the distributed benchmark. Neither this dataset DOI nor the software DOI is
a journal-article DOI; a future paper remains a separate research output.

## Official guidance

- [Enable GitHub archiving](https://help.zenodo.org/docs/github/enable-repository/)
- [CITATION.cff metadata](https://help.zenodo.org/docs/github/describe-software/citation-file/)
- [Archive a GitHub release](https://help.zenodo.org/docs/github/archive-software/github-upload/)
- [Reserve a dataset DOI in a draft](https://help.zenodo.org/docs/deposit/describe-records/reserve-doi/)

