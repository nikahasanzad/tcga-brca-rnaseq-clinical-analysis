# Final manifest policy

SHA256SUMS covers every public file except itself. FILE_INVENTORY.csv and
SOURCE_PROVENANCE_MANIFEST.csv cover all files except SHA256SUMS, themselves,
LICENSE_SCOPE.csv and RC2_TO_RC3_FILE_CHANGES.csv (self-referential ledgers).
LICENSE_SCOPE.csv covers every public path. RC2_TO_RC3_FILE_CHANGES.csv lists
non-ledger changes and retired historical reports. Root SHA256SUMS covers the
complete private/public local package except itself. Private source preservation
and predecessor manifests preserve the full lineage. RC2 is the exact audited
predecessor; unchanged records trace to its source-to-v09/V4 manifests.
New closure artifacts trace to private resource hashes, tests and the assembly script.
