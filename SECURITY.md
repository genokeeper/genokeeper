# Security

Genokeeper runs on your own machine against your own genome, so a security problem here is
usually a privacy problem: a way for genome data to leave the machine, to reach git, or for
code from a pull request to run where genome data is.

## Reporting

Report a vulnerability privately through GitHub's "Report a vulnerability" button on the
repository's Security tab, not in a public issue. Never attach genome data, findings or
outputs of `genokeeper` to a report: describe the problem on a public sample (GIAB's HG002) or
on invented data.

## Scope

- anything that sends reads, genotypes or findings off the machine without the user typing a
  command that says so (README, "Privacy and network");
- a way past the genome-data guard (`scripts/check_no_genome_data.sh`);
- a download the bootstrap does not check by checksum, or a container that runs with a
  network;
- the CI workflows (`.github/workflows/`).

## For maintainers

Review a pull request from outside in a checkout with no genome data and no lab: the git hooks,
`uv run genokeeper` and an agent's project settings all run the checked-out code.
