# Execution and offline checks

From repository/, run the self-contained base-R suite:

~~~sh
Rscript --vanilla tools/run_offline_tests.R .
~~~

It verifies SHA256, frozen syntax, synthetic runtime/QC, archived table counts and
inert help. No TCGA model is executed. The text fixture is public and self-contained.
With recorded graphics dependencies, also run:

~~~sh
Rscript --vanilla tests/test_plot_snapshot.R .
Rscript --vanilla tools/check_environment.R
~~~

For a separately authorized real-data run, arrange sibling private_inputs,
private_cache and private_runs directories outside this repository. All must
exist; input and output roots must be separate and non-nested. Use a new run ID.
The runner accepts relative input/run paths but requires an absolute cache path:

~~~sh
Rscript --vanilla run.R \
  --input-dir ../../private_inputs \
  --cache-dir "$(cd ../../private_cache && pwd -P)" \
  --runs-dir ../../private_runs \
  --run-id rc2-authorized-new-run \
  --execute
~~~

This command was not executed during closure. Offline is default. Do not add the
network flag when reproducing archived inputs; a current live acquisition has
separate provenance. The CLI refuses an existing run destination.
See INPUT_DATA.md and environment/README.md for the exact input/environment gates.
