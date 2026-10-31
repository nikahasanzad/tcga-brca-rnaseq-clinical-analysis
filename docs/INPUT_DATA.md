# Data acquisition and historical inputs

The source is NCI GDC TCGA-BRCA, STAR - Counts, Primary Tumor gene-expression
quantification. See [GDC API](https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/),
[file downloads](https://docs.gdc.cancer.gov/API/Users_Guide/Downloading_Files/) and
[Data Transfer Tool](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/).

data-manifests/gdc_final_manifest.tsv contains the exact 1,093 selected expression
file UUIDs and archived MD5/size/state. gdc_technical_qc_manifest.tsv contains the
ten paired-aliquot QC files; five also occur in the selected final set.
All originate from archived access=open metadata. No data download was performed
in this release closure, and current availability is not guaranteed.

For a separately authorized download, using an installed GDC client:

~~~sh
gdc-client download -m data-manifests/gdc_final_manifest.tsv -d ../../private_inputs/data_raw/rna_seq_final
gdc-client download -m data-manifests/gdc_technical_qc_manifest.tsv -d ../../private_inputs/data_raw/rna_seq_technical_qc
~~~

Preserve original file UUIDs and MD5; do not request latest replacements and call
them the same historical snapshot. Verify downloaded bytes before execution.
Controlled access, if encountered, needs separate authorization; no token is bundled.

The runner also requires data_raw/clinical/clinical.tsv and cache envelopes
query_exp.rds and biospecimen_raw.rds. Original API responses were unavailable.
For exact historical caches, obtain the six archived metadata exports listed in
data-manifests/metadata_input_fingerprints.csv from the authorized custodian.
That list also fingerprints the clinical TSV. They are intentionally not distributed.
The August 2026 export is distinct from a new current GDC export.

Rebuild metadata caches into a new directory, without reading analytical results:

~~~sh
Rscript --vanilla tools/prepare_archived_cache.R --input-dir ../../private_export_snapshot --cache-dir ../../private_cache_rc2 --execute
~~~

The portable builder requires jsonlite 2.0.0 and refuses an existing destination.
Its archived join logic is preserved. Closure tested the rebuilt envelopes against
the archived cache objects; see the metadata-cache QA record.
This is metadata reconstruction, not RNA-seq reanalysis or a new live-GDC test.

Full-data reproduction therefore still needs authorized access to the historical
clinical/export snapshot. Offline repository tests need none of these private files.
