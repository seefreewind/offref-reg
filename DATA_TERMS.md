# Data and output terms

**AlphaGenome outputs and derivatives in `source_data/`, `qa/frozen_bundle.json`, `tables/`, `supplementary/` and `figures/` are subject to the [AlphaGenome Output Terms of Use](https://deepmind.google.com/science/alphagenome/output-terms). By using this information, you agree to those terms. They are not licensed under MIT.**

The [services terms](https://deepmind.google.com/science/alphagenome/terms) permit publication/sharing with conspicuous notice. Output use is restricted to non-commercial research and remains subject to the applicable restrictions on model training, clinical use and misleading representation. Consult the linked binding terms; the terms effective on the output-generation date apply. This audit's Phase1R3 outputs were generated in October 2026; output terms were verified on 14 October 2026 (effective 9 October 2026).

Modifications: paired REF/ALT predictions were integrated over homologous target-gene exon intersections, differenced, harmonized, summarized and visualized. The package includes scores, statistical summaries and figures, not raw genome-length prediction arrays. Cross-model agreement is not biological validation.

GEUVADIS reference beta/PIP subsets derive from eQTL Catalogue release 7, QTD000110, https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000013/QTD000110/QTD000110.credible_sets.tsv.gz. These data are CC BY 4.0 under https://www.ebi.ac.uk/eqtl/License/. Attribute eQTL Catalogue (Kerimov et al., Nature Genetics 2021; PLoS Genetics 2023) and the GEUVADIS study. Our modifications were outcome-independent sampling, target mapping and prediction-score joining.

Insertion identities, aggregate annotations and project-derived reliability summaries originate from the HGSVC3 2024 v1.0 discovery panel and GNRS accessions PRJCA007976/GVM000324/GWHBHSK00000000. No individual genotypes, donor identifiers, expression matrices, third-party FASTA/VCF files, downloaded annotation databases, model weights or credentials are redistributed. Obtain original inputs from their providers under their own terms. This release does not grant rights to those upstream resources.

Borzoi predictions were generated using the official model ecosystem (https://github.com/calico/borzoi). Its code license must not be interpreted as a new license for third-party weights. No upstream code/weights are vendored here.

All data are available for verification of this non-commercial research subject to the above notices. No blanket open-data license overrides upstream restrictions. Original code is MIT.
