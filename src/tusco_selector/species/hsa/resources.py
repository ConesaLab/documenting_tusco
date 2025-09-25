"""
Human (Homo sapiens) specific resources.
"""

EXTERNAL_RESOURCES = {
    "annotation": {
        "gencode": {
            "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz",
            "output": "gencode.v49.annotation.gtf.gz",
            "description": "Gencode v49 annotation",
        },
        "refseq": {
            "refseq_annotation": {
                "url": "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/current/GCF_000001405.40-RS_2025_08/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz",
                "output": "GCF_000001405.40_GRCh38.p14_genomic.gtf.gz",
                "description": "RefSeq Annotation GCF_000001405.40",
            },
            "refseq_report": {
                "url": "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/current/GCF_000001405.40-RS_2025_08/GCF_000001405.40_GRCh38.p14_assembly_report.txt",
                "output": "GCF_000001405.40_GRCh38.p14_assembly_report.txt",
                "description": "RefSeq GCF_000001405.40_GRCh38.p14 assembly report",
            },
        },
        "mane": {
            "url": "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz",
            "output": "mane.v1.4.ensembl_genomic.gtf.gz",
            "description": "MANE v1.4 Ensembl Genomic",
        },
        "mapping": {
            "ensembl_refseq_biomart": {
                "url": (
                    "https://www.ensembl.org/biomart/martservice?"
                    "query=%3C%21DOCTYPE%20Query%3E%0A"
                    "%3CQuery%20virtualSchemaName%3D%22default%22%20formatter%3D%22TSV%22%20header%3D%220%22%0A"
                    "%20%20%20%20%20%20%20%20uniqueRows%3D%220%22%20count%3D%22%22%20datasetConfigVersion%3D%220.6%22%3E%0A"
                    "%20%20%3CDataset%20name%3D%22hsapiens_gene_ensembl%22%20interface%3D%22default%22%3E%0A"
                    "%20%20%20%20%3CAttribute%20name%3D%22ensembl_gene_id%22%2F%3E%0A"
                    "%20%20%20%20%3CAttribute%20name%3D%22ensembl_transcript_id%22%2F%3E%0A"
                    "%20%20%20%20%3CAttribute%20name%3D%22external_gene_name%22%2F%3E%0A"
                    "%20%20%20%20%3CAttribute%20name%3D%22entrezgene_id%22%2F%3E%0A"
                    "%20%20%20%20%3CAttribute%20name%3D%22refseq_mrna%22%2F%3E%0A"
                    "%20%20%20%20%3CAttribute%20name%3D%22refseq_peptide%22%2F%3E%0A"
                    "%20%20%3C%2FDataset%3E%0A"
                    "%3C%2FQuery%3E"
                ),
                "output": "ensembl_refseq.tsv",
                "description": "Ensembl gene + transcript ↔ RefSeq ↔ NCBI GeneID (human, GRCh38)",
            }
        },
    },
    "expression": {
        "gtex": {
            "adipose_subcutaneous": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_adipose_subcutaneous.gct.gz",
                "output": "UBERON:0002190.gct.gz",
                "description": "GTEx v10 Adipose - Subcutaneous",
            },
            "adipose_visceral_omentum": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_adipose_visceral_omentum.gct.gz",
                "output": "UBERON:0010414.gct.gz",
                "description": "GTEx v10 Adipose - Visceral (Omentum )",
            },
            "adrenal_gland": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_adrenal_gland.gct.gz",
                "output": "UBERON:0002369.gct.gz",
                "description": "GTEx v10 Adrenal Gland",
            },
            "artery_aorta": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_artery_aorta.gct.gz",
                "output": "UBERON:0001496.gct.gz",
                "description": "GTEx v10 Artery - Aorta",
            },
            "artery_coronary": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_artery_coronary.gct.gz",
                "output": "UBERON:0001621.gct.gz",
                "description": "GTEx v10 Artery - Coronary",
            },
            "artery_tibial": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_artery_tibial.gct.gz",
                "output": "UBERON:0007610.gct.gz",
                "description": "GTEx v10 Artery - Tibial",
            },
            "bladder": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_bladder.gct.gz",
                "output": "UBERON:0001255.gct.gz",
                "description": "GTEx v10 Bladder",
            },
            "brain_amygdala": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_brain_amygdala.gct.gz",
                "output": "UBERON:0001876.gct.gz",
                "description": "GTEx v10 Brain - Amygdala",
            },
            "brain_anterior_cingulate_cortex_ba24": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_brain_anterior_cingulate_cortex_ba24.gct.gz",
                "output": "UBERON:0009835.gct.gz",
                "description": "GTEx v10 Brain - Anterior cingulate cortex (BA24 )",
            },
            "brain_caudate_basal_ganglia": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_brain_caudate_basal_ganglia.gct.gz",
                "output": "UBERON:0001873.gct.gz",
                "description": "GTEx v10 Brain - Caudate (basal ganglia )",
            },
            "brain_cerebellar_hemisphere": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_brain_cerebellar_hemisphere.gct.gz",
                "output": "UBERON:0002037.gct.gz",
                "description": "GTEx v10 Brain - Cerebellar Hemisphere",
            },
            "brain_cerebellum": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_brain_cerebellum.gct.gz",
                "output": "UBERON:0002037.gct.gz",
                "description": "GTEx v10 Brain - Cerebellum",
            },
            "brain_cortex": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_brain_cortex.gct.gz",
                "output": "UBERON:0001870.gct.gz",
                "description": "GTEx v10 Brain - Cortex",
            },
            "brain_frontal_cortex_ba9": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_brain_frontal_cortex_ba9.gct.gz",
                "output": "UBERON:0009834.gct.gz",
                "description": "GTEx v10 Brain - Frontal Cortex (BA9 )",
            },
            "brain_hippocampus": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_brain_hippocampus.gct.gz",
                "output": "UBERON:0001954.gct.gz",
                "description": "GTEx v10 Brain - Hippocampus",
            },
            "brain_hypothalamus": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_brain_hypothalamus.gct.gz",
                "output": "UBERON:0001898.gct.gz",
                "description": "GTEx v10 Brain - Hypothalamus",
            },
            "brain_nucleus_accumbens_basal_ganglia": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_brain_nucleus_accumbens_basal_ganglia.gct.gz",
                "output": "UBERON:0001882.gct.gz",
                "description": "GTEx v10 Brain - Nucleus accumbens (basal ganglia )",
            },
            "brain_putamen_basal_ganglia": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_brain_putamen_basal_ganglia.gct.gz",
                "output": "UBERON:0001874.gct.gz",
                "description": "GTEx v10 Brain - Putamen (basal ganglia )",
            },
            "brain_spinal_cord_cervical_c-1": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_brain_spinal_cord_cervical_c-1.gct.gz",
                "output": "UBERON:0006469.gct.gz",
                "description": "GTEx v10 Brain - Spinal cord (cervical c-1 )",
            },
            "brain_substantia_nigra": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_brain_substantia_nigra.gct.gz",
                "output": "UBERON:0002038.gct.gz",
                "description": "GTEx v10 Brain - Substantia nigra",
            },
            "breast_mammary_tissue": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_breast_mammary_tissue.gct.gz",
                "output": "UBERON:0008367.gct.gz",
                "description": "GTEx v10 Breast - Mammary Tissue",
            },
            "cells_ebv-transformed_lymphocytes": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_cells_ebv-transformed_lymphocytes.gct.gz",
                "output": "EFO:0000572.gct.gz",
                "description": "GTEx v10 Cells - EBV-transformed lymphocytes",
            },
            "cells_cultured_fibroblasts": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_cells_cultured_fibroblasts.gct.gz",
                "output": "EFO:0002009.gct.gz",
                "description": "GTEx v10 Cells - Cultured fibroblasts",
            },
            "cervix_ectocervix": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_cervix_ectocervix.gct.gz",
                "output": "UBERON:0012249.gct.gz",
                "description": "GTEx v10 Cervix - Ectocervix",
            },
            "cervix_endocervix": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_cervix_endocervix.gct.gz",
                "output": "UBERON:0000458.gct.gz",
                "description": "GTEx v10 Cervix - Endocervix",
            },
            "colon_sigmoid": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_colon_sigmoid.gct.gz",
                "output": "UBERON:0001159.gct.gz",
                "description": "GTEx v10 Colon - Sigmoid",
            },
            "colon_transverse": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_colon_transverse.gct.gz",
                "output": "UBERON:0001157.gct.gz",
                "description": "GTEx v10 Colon - Transverse",
            },
            "esophagus_gastroesophageal_junction": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_esophagus_gastroesophageal_junction.gct.gz",
                "output": "UBERON:0004550.gct.gz",
                "description": "GTEx v10 Esophagus - Gastroesophageal Junction",
            },
            "esophagus_mucosa": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_esophagus_mucosa.gct.gz",
                "output": "UBERON:0006920.gct.gz",
                "description": "GTEx v10 Esophagus - Mucosa",
            },
            "esophagus_muscularis": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_esophagus_muscularis.gct.gz",
                "output": "UBERON:0004648.gct.gz",
                "description": "GTEx v10 Esophagus - Muscularis",
            },
            "fallopian_tube": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_fallopian_tube.gct.gz",
                "output": "UBERON:0003889.gct.gz",
                "description": "GTEx v10 Fallopian Tube",
            },
            "heart_atrial_appendage": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_heart_atrial_appendage.gct.gz",
                "output": "UBERON:0006631.gct.gz",
                "description": "GTEx v10 Heart - Atrial Appendage",
            },
            "heart_left_ventricle": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_heart_left_ventricle.gct.gz",
                "output": "UBERON:0006566.gct.gz",
                "description": "GTEx v10 Heart - Left Ventricle",
            },
            "kidney_cortex": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_kidney_cortex.gct.gz",
                "output": "UBERON:0001225.gct.gz",
                "description": "GTEx v10 Kidney - Cortex",
            },
            "kidney_medulla": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_kidney_medulla.gct.gz",
                "output": "UBERON:0001293.gct.gz",
                "description": "GTEx v10 Kidney - Medulla",
            },
            "liver": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_liver.gct.gz",
                "output": "UBERON:0001114.gct.gz",
                "description": "GTEx v10 Liver",
            },
            "lung": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_lung.gct.gz",
                "output": "UBERON:0008952.gct.gz",
                "description": "GTEx v10 Lung",
            },
            "minor_salivary_gland": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_minor_salivary_gland.gct.gz",
                "output": "UBERON:0006330.gct.gz",
                "description": "GTEx v10 Minor Salivary Gland",
            },
            "muscle_skeletal": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_muscle_skeletal.gct.gz",
                "output": "UBERON:0011907.gct.gz",
                "description": "GTEx v10 Muscle - Skeletal",
            },
            "nerve_tibial": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_nerve_tibial.gct.gz",
                "output": "UBERON:0001323.gct.gz",
                "description": "GTEx v10 Nerve - Tibial",
            },
            "ovary": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_ovary.gct.gz",
                "output": "UBERON:0000992.gct.gz",
                "description": "GTEx v10 Ovary",
            },
            "pancreas": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_pancreas.gct.gz",
                "output": "UBERON:0001150.gct.gz",
                "description": "GTEx v10 Pancreas",
            },
            "pituitary": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_pituitary.gct.gz",
                "output": "UBERON:0000007.gct.gz",
                "description": "GTEx v10 Pituitary",
            },
            "prostate": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_prostate.gct.gz",
                "output": "UBERON:0002367.gct.gz",
                "description": "GTEx v10 Prostate",
            },
            "skin_not_sun_exposed_suprapubic": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_skin_not_sun_exposed_suprapubic.gct.gz",
                "output": "UBERON:0036149.gct.gz",
                "description": "GTEx v10 Skin - Not Sun Exposed (Suprapubic )",
            },
            "skin_sun_exposed_lower_leg": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_skin_sun_exposed_lower_leg.gct.gz",
                "output": "UBERON:0004264.gct.gz",
                "description": "GTEx v10 Skin - Sun Exposed (Lower leg )",
            },
            "small_intestine_terminal_ileum": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_small_intestine_terminal_ileum.gct.gz",
                "output": "UBERON:0001211.gct.gz",
                "description": "GTEx v10 Small Intestine - Terminal Ileum",
            },
            "spleen": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_spleen.gct.gz",
                "output": "UBERON:0002106.gct.gz",
                "description": "GTEx v10 Spleen",
            },
            "stomach": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_stomach.gct.gz",
                "output": "UBERON:0000945.gct.gz",
                "description": "GTEx v10 Stomach",
            },
            "testis": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_testis.gct.gz",
                "output": "UBERON:0000473.gct.gz",
                "description": "GTEx v10 Testis",
            },
            "thyroid": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_thyroid.gct.gz",
                "output": "UBERON:0002046.gct.gz",
                "description": "GTEx v10 Thyroid",
            },
            "uterus": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_uterus.gct.gz",
                "output": "UBERON:0000995.gct.gz",
                "description": "GTEx v10 Uterus",
            },
            "vagina": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_vagina.gct.gz",
                "output": "UBERON:0000996.gct.gz",
                "description": "GTEx v10 Vagina",
            },
            "whole_blood": {
                "url": "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_whole_blood.gct.gz",
                "output": "UBERON:0013756.gct.gz",
                "description": "GTEx v10 Whole Blood",
            },
        },
        "hrt_atlas": {
            "url": "https://housekeeping.unicamp.br/Housekeeping_GenesHuman.csv",
            "output": "hrt_atlas.housekeeping_genes.csv",
            "description": "HRT Atlas Housekeeping Genes",
        },
        "bgee": {
            "bgee_expression": {
                "url": "https://www.bgee.org/ftp/current/download/calls/expr_calls/Homo_sapiens_expr_simple.tsv.gz",
                "output": "bgee.expr_simple.tsv.gz",
                "description": "Bgee human Expression",
            },
            "bgee_mapping": {
                "url": "https://www.bgee.org/ftp/current/download/processed_expr_values/rna_seq/Homo_sapiens/Homo_sapiens_RNA-Seq_experiments_libraries.tar.gz",
                "output": "bgee.rna_seq_experiments_libraries.tar.gz",
                "description": "Bgee human RNA-Seq Experiments and Libraries",
            },
        },
    },
    "splicing_tss": {
        "recount3": {
            "url": "https://ftp.ebi.ac.uk/pub/databases/havana/ngs_havana/human_jbrowse/gencode_primary/tracks/recount3/recount3_final.pass1V1.bed.gz",
            "output": "recount3.pass1V1.bed.gz",
            "description": "Recount3 Pass1 V1",
        },
        "refTSS": {
            "url": "http://reftss.clst.riken.jp/datafiles/current/human/refTSS_v4.1_human_coordinate.hg38.bed.txt.gz",
            "output": "refTSS.v4.1.human.hg38.bed.gz",
            "description": "RefTSS v4.1 Human HG38",
        },
        "introverse": {
            "url": "https://tusco.s3.us-east-1.amazonaws.com/IntroVerse_80.csv",
            "output": "introverse_80.csv",
            "description": "Filters genes present novel isoforms in >80% of samples",
        },
    },
    "manual_filter": {
        "url": "https://tusco.s3.us-east-1.amazonaws.com/manual_filtered.txt",
        "output": "manual_filtered.txt",
        "description": "Manual filtering of TUSCO dataset",
    },
}
