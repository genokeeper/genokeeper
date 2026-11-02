# Runbook: bug reports

How to file, close and reference an entry in this project's ledger. The index links
here; the entry files are the ledger itself.

File with `stayfixed bugs new`; close an entry by setting its status and `fixed_in`, then run
`stayfixed bugs index` rather than editing the index by hand. `stayfixed bugs check` is the gate
over every rule the ledger holds, and `stayfixed bugs renumber` moves an entry and rewrites
every mention of it. The `file-bug` and `close-bug` skills walk both ends of that.

**There is no merge.** Two entries that turn out to be one bug stay two files: close the later
one and name the survivor in its `related` list, which `stayfixed bugs check` holds to an entry
that exists. Nothing folds one entry into another — `bugs renumber` moves an entry to a free
identifier and rewrites the mentions of it, which is a different operation.
