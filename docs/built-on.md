# Built on

Genokeeper runs these tools and reads these truth sets and knowledge sources, in pipeline
order. The papers describe the tools, not Genokeeper's own accuracy, and no group listed
here endorses its output. What Genokeeper measured itself, and on which samples, is in
[reference.md](reference.md) and on the HG002 demo's
[calibration page](https://genokeeper.github.io/genokeeper/calibration.html).

**Reads to calls.**
- [GRCh38 reference assembly](https://doi.org/10.1101/gr.213611.116): Schneider et al., Genome Res 2017.
- [GIAB GRCh38 with GRC false duplications masked](https://doi.org/10.1038/s41587-021-01158-1): Wagner et al., Nat Biotechnol 2022.
- [BWA-MEM](https://doi.org/10.48550/arXiv.1303.3997): Li, arXiv 2013.
- [SAMtools](https://doi.org/10.1093/gigascience/giab008): Danecek et al., GigaScience 2021.
- [BCFtools](https://doi.org/10.1093/gigascience/giab008): Danecek et al., GigaScience 2021.
- [mosdepth](https://doi.org/10.1093/bioinformatics/btx699): Pedersen and Quinlan, Bioinformatics 2018.
- [VerifyBamID2](https://doi.org/10.1101/gr.246934.118): Zhang et al., Genome Res 2020.
- [DeepVariant](https://doi.org/10.1038/nbt.4235): Poplin et al., Nat Biotechnol 2018.
- [mutserve (mtDNA-Server)](https://doi.org/10.1093/nar/gkw247): Weissensteiner et al., Nucleic Acids Res 2016.
- [WhatsHap](https://doi.org/10.1101/085050): Martin et al., bioRxiv 2016.

**Truth sets and calibration.**
- [Genome in a Bottle reference samples (HG002, HG003, HG004)](https://doi.org/10.1038/sdata.2016.25): Zook et al., Sci Data 2016.
- [GIAB v4.2.1 small-variant benchmarks (HG003, HG004)](https://doi.org/10.1016/j.xgen.2022.100128): Wagner et al., Cell Genom 2022.
- [GIAB HG002 v5.0q benchmark (T2T Q100 v1.1 assembly)](https://doi.org/10.1101/2025.09.21.677443): Hansen et al., bioRxiv 2025.
- [GIAB genome stratifications v3.6](https://doi.org/10.1038/s41467-024-53260-y): Dwarshuis et al., Nat Commun 2024.
- [RTG Tools vcfeval](https://doi.org/10.1101/023754): Cleary et al., bioRxiv 2015.

**Structural variants, copy number and repeats.**
- [DELLY](https://doi.org/10.1093/bioinformatics/bts378): Rausch et al., Bioinformatics 2012.
- [Dysgu](https://doi.org/10.1093/nar/gkac039): Cleal and Baird, Nucleic Acids Res 2022.
- [Truvari](https://doi.org/10.1186/s13059-022-02840-6): English et al., Genome Biol 2022.
- [CNVpytor](https://doi.org/10.1093/gigascience/giab074): Suvakov et al., GigaScience 2021.
- [INSurVeyor](https://doi.org/10.1038/s41467-023-38870-2): Rajaby et al., Nat Commun 2023.
- [SCRAMble](https://doi.org/10.1038/s41436-020-0749-x): Torene et al., Genet Med 2020.
- [ExpansionHunter](https://doi.org/10.1093/bioinformatics/btz431): Dolzhenko et al., Bioinformatics 2019.
- [REViewer](https://doi.org/10.1186/s13073-022-01085-z): Dolzhenko et al., Genome Med 2022.
- [STRchive](https://doi.org/10.1186/s13073-025-01454-4): Hiatt et al., Genome Med 2025.
- [SMNCopyNumberCaller](https://doi.org/10.1038/s41436-020-0754-0): Chen et al., Genet Med 2020.
- [Parascopy](https://doi.org/10.1038/s41467-022-30930-3): Prodanov and Bansal, Nat Commun 2022.
- [Gauchian](https://doi.org/10.1038/s42003-022-03610-7): Toffoli et al., Commun Biol 2022.
- [KILDA](https://doi.org/10.1093/nargab/lqaf070): Molitor et al., NAR Genom Bioinform 2025.

**Pharmacogenomics.**
- [PharmCAT](https://doi.org/10.1002/cpt.1568): Sangkuhl et al., Clin Pharmacol Ther 2020.
- [Aldy](https://doi.org/10.1101/gr.277075.122): Hari et al., Genome Res 2023.
- [BCyrius](https://doi.org/10.1002/prp2.70065): Halman and Conyers, Pharmacol Res Perspect 2025.
- [StellarPGx](https://doi.org/10.1002/cpt.2173): Twesigomwe et al., Clin Pharmacol Ther 2021.
- [CPIC](https://doi.org/10.1038/clpt.2010.279): Relling and Klein, Clin Pharmacol Ther 2011.
- [PharmGKB / ClinPGx clinical annotations](https://doi.org/10.1002/cpt.2350): Whirl-Carrillo et al., Clin Pharmacol Ther 2021.

**Clinical and population evidence.**
- [ClinVar](https://doi.org/10.1093/nar/gkae1090): Landrum et al., Nucleic Acids Res 2025.
- [ClinGen](https://doi.org/10.1056/NEJMsr1406261): Rehm et al., N Engl J Med 2015.
- [ClinGen gene-disease validity](https://doi.org/10.1016/j.ajhg.2017.04.015): Strande et al., Am J Hum Genet 2017.
- [ClinGen dosage sensitivity](https://search.clinicalgenome.org/kb/gene-dosage): ClinGen Dosage Sensitivity Curation.
- [ClinGen SVI BA1 exception list](https://doi.org/10.1002/humu.23642): Ghosh et al., Hum Mutat 2018.
- [GenCC](https://doi.org/10.1016/j.gim.2022.04.017): DiStefano et al., Genet Med 2022.
- [Clinical Genomic Database (CGD)](https://doi.org/10.1073/pnas.1302575110): Solomon et al., Proc Natl Acad Sci USA 2013.
- [ACMG SF v3.3 secondary findings list](https://doi.org/10.1016/j.gim.2025.101454): Lee et al., Genet Med 2025.
- [ACMG carrier screening tiers (2021)](https://doi.org/10.1038/s41436-021-01203-z): Gregg et al., Genet Med 2021.
- [gnomAD v4.1 (allele frequencies and SV sites)](https://doi.org/10.1038/s41586-023-06045-0): Chen et al., Nature 2024.
- [Filtering allele frequency (FAF)](https://doi.org/10.1038/gim.2017.26): Whiffin et al., Genet Med 2017.
- [MANE transcripts](https://doi.org/10.1038/s41586-022-04558-8): Morales et al., Nature 2022.
- [MITOMAP](https://doi.org/10.1002/0471250953.bi0123s44): Lott et al., Curr Protoc Bioinformatics 2013.
- [Genomics England PanelApp](https://doi.org/10.1038/s41588-019-0528-2): Martin et al., Nat Genet 2019.
- [GWAS Catalog](https://doi.org/10.1093/nar/gkae1070): Cerezo et al., Nucleic Acids Res 2025.
- [LitVar 2.0](https://doi.org/10.1038/s41588-023-01414-x): Allot et al., Nat Genet 2023.
- [MedlinePlus Genetics](https://medlineplus.gov/genetics/): National Library of Medicine, MedlinePlus Genetics.

**Variant effect prediction.**
- [Ensembl VEP](https://doi.org/10.1186/s13059-016-0974-4): McLaren et al., Genome Biol 2016.
- [SpliceAI](https://doi.org/10.1016/j.cell.2018.12.015): Jaganathan et al., Cell 2019.
- [REVEL](https://doi.org/10.1016/j.ajhg.2016.08.016): Ioannidis et al., Am J Hum Genet 2016.
- [AlphaMissense](https://doi.org/10.1126/science.adg7492): Cheng et al., Science 2023.
- [SpliceVault](https://doi.org/10.1038/s41588-022-01293-8): Dawes et al., Nat Genet 2023.
- [pext (proportion expressed across transcripts)](https://doi.org/10.1038/s41586-020-2329-2): Cummings et al., Nature 2020.
- [UCSC NMD Escape track (MANE)](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=nmd): UCSC Genome Browser, NMD Escape track.

**Polygenic scores and ancestry.**
- [PGS Catalog](https://doi.org/10.1038/s41588-021-00783-5): Lambert et al., Nat Genet 2021.
- [pgsc_calc](https://doi.org/10.1038/s41588-024-01937-x): Lambert et al., Nat Genet 2024.
- [Nextflow](https://doi.org/10.1038/nbt.3820): Di Tommaso et al., Nat Biotechnol 2017.
- [PLINK 2](https://doi.org/10.1186/s13742-015-0047-8): Chang et al., GigaScience 2015.
- [FRAPOSA](https://doi.org/10.1093/bioinformatics/btaa152): Zhang et al., Bioinformatics 2020.
- [HGDP + 1000 Genomes reference panel (gnomAD v3.1)](https://doi.org/10.1101/gr.278378.123): Koenig et al., Genome Res 2024.
- [fastmixture](https://doi.org/10.24072/pcjournal.503): Santander et al., Peer Community J 2024.
- [hmmix](https://doi.org/10.1371/journal.pgen.1007641): Skov et al., PLoS Genet 2018.

**Haplogroups, HLA and other typing.**
- [Haplogrep 3](https://doi.org/10.1093/nar/gkad284): Schönherr et al., Nucleic Acids Res 2023.
- [PhyloTree Build 17](https://doi.org/10.1016/j.fsigss.2015.09.155): van Oven, Forensic Sci Int Genet Suppl Ser 2015.
- [haplocheck](https://doi.org/10.1101/gr.256545.119): Weissensteiner et al., Genome Res 2021.
- [Yleaf](https://doi.org/10.1093/molbev/msy032): Ralf et al., Mol Biol Evol 2018.
- [OptiType](https://doi.org/10.1093/bioinformatics/btu548): Szolek et al., Bioinformatics 2014.
- [T1K](https://doi.org/10.1101/gr.277585.122): Song et al., Genome Res 2023.
- [IPD-IMGT/HLA](https://doi.org/10.1093/nar/gkac1011): Barker et al., Nucleic Acids Res 2023.
- [RBCeq2](https://github.com/limcintyre/RBCeq2): McIntyre et al., RBCeq2 v2.4.4 (software) 2026.
- [BCFtools/RoH](https://doi.org/10.1093/bioinformatics/btw044): Narasimhan et al., Bioinformatics 2016.
