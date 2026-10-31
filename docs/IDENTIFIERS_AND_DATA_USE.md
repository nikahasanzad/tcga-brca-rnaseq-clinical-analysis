# TCGA identifiers and source data conditions

The fixed code labels are TCGA participant, specimen or aliquot submitter identifiers;
they are not names, email addresses or a key to personal identities. The archived
repository metadata contains 1,111 expression records marked access=open and
state=released. The fixed labels are used to reproduce specimen/aliquot selection,
diagnostic inspection and six historical exclusions. They are retained unchanged.

Based on their open metadata origin and provenance function, retaining these labels
in reproducible code is appropriate for this candidate's traceability purpose.
This is a scoped release assessment, not blanket authorization to distribute all
local patient-level derivatives. rc1's categorical identifier blocker is replaced
by this source-specific assessment. No attempt to identify participants was made.

[GDC policy](https://gdc.cancer.gov/access-data/data-access-policies) distinguishes
open from controlled data, requires acknowledgment and prohibits re-identification.
[The API guide](https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/)
explains that file/case/sample UUIDs are resource identifiers.
Public download manifests contain file UUIDs, filenames, MD5, size and archived
state only; patient mappings and phenotypes are omitted.

Public policy review does not establish a current live access check for every
historical file. Readers must check current access and source conditions before
downloading. The owner approved scoped project licensing; eventual online publication remains a separate action.
Raw clinical tables, detailed deletion records, raw logs and caches remain local.
