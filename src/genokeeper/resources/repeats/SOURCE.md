# The repeat catalog cut

`stage3a.GRCh38.json` is the variant catalog ExpansionHunter and REViewer read for the repeat
loci (the `locus_id`s of `[repeat_loci]` in config/grading.toml: slice 3a's eleven and slice
3e's thirteen) and for RFC1 and CFTR. The name is slice 3a's, kept so the paths that name it
stay valid.

- Source: Illumina's ExpansionHunter GRCh38 variant catalog at tag v5.0.0 (the version the
  image runs), `variant_catalog/grch38/variant_catalog.json`,
  https://github.com/Illumina/ExpansionHunter/blob/v5.0.0/variant_catalog/grch38/variant_catalog.json
- sha256 of the source file:
  `41043d6cac8a535922c3e29960991e09aaad3b24b14595a7d864232c455781a9`
- Why Illumina's and not STRchive's (design decision, 2026-10-28): STRchive v2.28.0's
  ExpansionHunter asset uses its locus bounds, which carry part of a unit at six loci and stop
  two CAG short at HTT; on them ExpansionHunter called HG002 one repeat above GIAB's v5.0q
  truth at DMPK, CACNA1A, ATXN1, ATXN2, ATXN7 and ATN1. On Illumina's regions it matched the
  truth at eight of the eleven loci. The thresholds stay STRchive's.
- The added entries, `added.GRCh38.json` (written by hand, `VariantType` `Repeat`): Illumina's
  catalog lacks ZFHX3, THAP11 and FGF14 and counts PABPN1's pure GCG run (`(GCG)*`, 18 bp, 6
  units), so these four come from STRchive v2.28.0's ExpansionHunter asset
  (`STRchive-disease-loci-v2.28.0.hg38.expansionhunter.json`, sha256
  `e87177dcae3e09acaa3fe766818e16d334b1f206f29f714bcf6b42fd2f5a5c9a`), each trimmed to whole
  units in its motif's phase (STRchive's carry 10.33, 21.33, 29.33 and 50.67 units, the flaw
  the 2026-10-28 decision rejected); CFTR's TG/T tract is read off GRCh38 and ClinVar
  2026-10-11 (VCV 229652, 496577, 818109, 1167538: `NC_000007.14:g.117548607_117548628` is
  TG×11 and `…_117548635` adds T×7; MANE NM_000492.4 c.1210-34TG[11]T[7], legacy IVS8).
  0-based start, end exclusive, as the catalog writes them:

  | Locus | Structure | Region | STRchive's region | Reference tract |
  |---|---|---|---|---|
  | PABPN1 | `(GCN)*` | `chr14:23321472-23321502` | `chr14:23321472-23321503` | GCG×6 GCA×3 GCG: 10 GCN, STRchive's frame |
  | ZFHX3 | `(CCG)*` | `chr16:72787695-72787755` | `chr16:72787694-72787758` | 20 units, CCA/CTG interruptions |
  | THAP11 | `(CAG)*` | `chr16:67842863-67842950` | `chr16:67842862-67842950` | 29 units, CAA interruptions |
  | FGF14 | `(AAG)*` | `chr13:102161574-102161724` | `chr13:102161574-102161726` | 50 AAG (= GAA on the gene) |
  | CFTR | `(TG)*(T)*` | `chr7:117548606-117548628`, `chr7:117548628-117548635` | — | TG×11 then T×7; VariantIds `CFTR_TG`, `CFTR_T` |

- AFF2 keeps Illumina's 20-unit region (`chrX:148500631-148500691`): STRchive's
  (`chrX:148500604-148500753`) reads 50 units in the reference, above STRchive's own benign
  maximum of 39.
- RFC1 is Illumina's entry (`(AARRG)*`, `chr4:39348424-39348479`) with
  `"RFC1MotifAnalysis": true`, which turns on the motif analysis of the ExpansionHunter build
  at commit 9e93322f4a; its grade is the analysis's call, not a `[repeat_loci]` row.
- The cut's rule (`cut_catalog`): Illumina's entries whose `LocusId` is a locus's, in the
  source's file order, each unchanged but for its `OfftargetRegions` (removed) and its
  `ReferenceRegion`s' contigs, named as the reference names them (`chr` prefixed), PABPN1's
  (the one `REPLACED` name) swapped in place for the added file's; then the added file's other
  entries of the loci, in that file's order; then RFC1's entry, then CFTR's. A `LocusId` in
  both files but not `REPLACED`, a locus in neither, and a missing RFC1 or CFTR entry are
  errors. The file is the rule's output as `json.dumps(entries, indent=2)` and a final
  newline.
- Regenerate: `uv run python -m genokeeper.repeats.expansionhunter cut <source file>`.

## Out

STRchive v2.28.0 lists 79 genes. Twenty-four are `[repeat_loci]` rows, RFC1 is graded by its
motif analysis, and 22 are written as not covered (`genokeeper.repeats.not_covered`). These 32
are left out, with the reason:

| Gene | Reason |
|---|---|
| AFF3 | congenital or childhood onset, no adult screening use |
| ARX | congenital or childhood onset, no adult screening use |
| CBL | congenital or childhood onset, no adult screening use |
| CEL | congenital or childhood onset, no adult screening use |
| COMP | congenital or childhood onset, no adult screening use |
| CSNK1E | congenital or childhood onset, no adult screening use |
| EIF4A3 | congenital or childhood onset, no adult screening use |
| FOXL2 | congenital or childhood onset, no adult screening use |
| GLS | congenital or childhood onset, no adult screening use |
| HOXA13 | congenital or childhood onset, no adult screening use |
| HOXD13 | congenital or childhood onset, no adult screening use |
| MIR7-2 | congenital or childhood onset, no adult screening use |
| NAXE | congenital or childhood onset, no adult screening use |
| PHOX2B | congenital or childhood onset, no adult screening use |
| PRDM12 | congenital or childhood onset, no adult screening use |
| RUNX2 | congenital or childhood onset, no adult screening use |
| SOX3 | congenital or childhood onset, no adult screening use |
| TBX1 | congenital or childhood onset, no adult screening use |
| TYMS | congenital or childhood onset, no adult screening use |
| VWA1 | congenital or childhood onset, no adult screening use |
| XYLT1 | congenital or childhood onset, no adult screening use |
| ZIC2 | congenital or childhood onset, no adult screening use |
| ZIC3 | congenital or childhood onset, no adult screening use |
| ZNF713 | congenital or childhood onset, no adult screening use |
| DIP2B | evidence disputed |
| DMD | evidence refuted |
| NIPA1 | evidence disputed |
| POLG | evidence disputed |
| MUC1 | a VNTR beyond short reads |
| PLIN4 | a VNTR beyond short reads |
| PRNP | an octarepeat beyond short reads |
| TAF1 | a retrotransposon-linked Filipino founder allele |
