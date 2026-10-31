# Self contained release tests

tools/run_offline_tests.R runs all data-free controls under base R and returns
nonzero on a failure. test_structure.R compares parsed frozen expressions to
tests/fixtures/scientific_contract.R, a text conversion proven identical to the
original private RDS syntax fixture. Neither syntax fixture evaluates the TCGA model.

test_runtime.R uses disposable synthetic data and mocked network; test_optional_qc.R
exercises a missing-value helper; test_public_tables.R checks saved result counts.
SHA256 validation covers the release contents. A negative tampering test in closure
confirms changed bytes fail the integrity gate.

test_plot_snapshot.R uses only synthetic graphics and recorded ggplot2.
The default GitLab job uses verified rocker/r-ver:4.5.1 by digest, needs no data,
and does not install scientific packages. Remote execution is a separate untested
step; local macOS success is not represented as a Linux runner result.
