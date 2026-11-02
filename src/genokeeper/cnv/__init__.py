"""Slice 3d: copy-number changes of 1 kb or more from read depth, in the screen genes.

`callers` runs the two depth callers on the CRAM (DELLY's `cnv` mode and CNVpytor) and reads
their losses and gains, and Dysgu's duplication junctions, as `deletions.screen.Deletion`
events of type `DEL` or `DUP`, which 3b's screen machinery takes.
"""
