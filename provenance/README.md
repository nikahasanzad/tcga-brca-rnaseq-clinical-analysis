# Source provenance and integrity

[SOURCE_PROVENANCE_MANIFEST.csv](../SOURCE_PROVENANCE_MANIFEST.csv) records
source aliases, exact source/current SHA256 and packaging transformations.
RC1 is the immediate copy lineage; analytical code additionally resolves directly
to audited CODE_V09. Original rc1 lineage/manifests remain in local_only/rc1_evidence.
EVIDENCE_INDEX.json and docs/EVIDENCE_INDEX.md resolve E01-E10 to actual bundled files.

Generated final ledgers/inventories exclude self-hashes to avoid recursion.
SHA256SUMS covers every other public repository file. FILE_INVENTORY.csv excludes
itself and the checksum file. Exact absolute source roots and preservation
inventories are private, never required for the public offline tests.
Hashes establish file identity, not scientific independence or cryptographic signatures.
