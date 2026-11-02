"""Slice 3b: deletions of 50 bp or more in the ACMG SF and carrier genes.

`files` names the inputs the bootstrap fetches (gnomAD-SV and HG002's v5.0q structural-variant
truth); `exons` cuts the MANE exons of the listed genes that the screen reads; `callers` runs
DELLY and Dysgu on the CRAM and reads their deletions; `truvari` matches the two callers'
deletions and, for the calibration, the calls against HG002's truth; `screen` shows or drops
each called deletion and puts it in its call class.
"""
