"""
Mouse (Mus musculus) specific resources.
"""

EXTERNAL_RESOURCES = {
    "annotation": {
        "gencode": {
            "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.annotation.gtf.gz",
            "output": "gencode.vM38.annotation.gtf.gz",
            "description": "Gencode vM38 annotation",
        },
        "refseq": {
            "refseq_annotation": {
                "url": "https://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/annotation_releases/current/GCF_000001635.27-RS_2024_02/GCF_000001635.27_GRCm39_genomic.gtf.gz",
                "output": "GCF_000001635.27_GRCm39_genomic.gtf.gz",
                "description": "RefSeq Annotation Release 100090",
            },
            "refseq_report": {
                "url": "https://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/annotation_releases/current/GCF_000001635.27-RS_2024_02/GCF_000001635.27_GRCm39_assembly_report.txt",
                "output": "GCF_000001635.27_GRCm39_assembly_report.txt",
                "description": "RefSeq Annotation Release 100090 assembly report",
            },
        },
        "reftss": {
            "url": "https://reftss.riken.jp/datafiles/current/mouse/refTSS_v4.1_mouse_coordinate.mm39.bed.gz",
            "output": "refTSS_v4.1_mouse_coordinate.mm39.bed.gz",
            "description": "Reference TSS v4.1 mouse coordinates (mm39)",
        },
        "recount3": {
            "url": "https://ftp.ebi.ac.uk/pub/databases/havana/recount3/m39_recount3_final_sorted.pass1.bed.gz",
            "output": "m39_recount3_final_sorted.pass1.bed.gz",
            "description": "Recount3 Pass1 Mouse",
        },
        "mapping": {
            "ensembl_refseq_biomart": {
                "url": (
                    "https://www.ensembl.org/biomart/martservice?"
                    "query=%3C%21DOCTYPE%20Query%3E%0A"
                    "%3CQuery%20virtualSchemaName%3D%22default%22%20formatter%3D%22TSV%22%20header%3D%220%22%0A"
                    "%20%20%20%20%20%20%20%20uniqueRows%3D%220%22%20count%3D%22%22%20datasetConfigVersion%3D%220.6%22%3E%0A"
                    "%20%20%3CDataset%20name%3D%22mmusculus_gene_ensembl%22%20interface%3D%22default%22%3E%0A"
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
                "description": "Ensembl gene + transcript ↔ RefSeq ↔ NCBI GeneID (mouse, GRCm39)",
            }
        },
        "housekeeping": {
            "housekeeping_genes": {
                "url": "https://housekeeping.unicamp.br/Housekeeping_GenesMouse.csv",
                "output": "Housekeeping_GenesMouse.csv",
                "description": "Housekeeping genes (mouse)",
            },
            "manual_filtered": {
                "url": "https://tusco.s3.us-east-1.amazonaws.com/manual_filtered_mouse.txt",
                "output": "manual_filtered_mouse.txt",
                "description": "Manually filtered mouse housekeeping genes",
            },
        },
    },
    "expression": {
        "encode": {
            "macrophage": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=macrophage",
                "output": "encode.macrophage.expression.tsv",
                "description": "Encode RNA expression - macrophage",
            },
            "heart": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=heart",
                "output": "encode.heart.expression.tsv",
                "description": "Encode RNA expression - heart",
            },
            "adrenal_gland": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=adrenal+gland",
                "output": "encode.adrenal_gland.expression.tsv",
                "description": "Encode RNA expression - adrenal gland",
            },
            "liver": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=liver",
                "output": "encode.liver.expression.tsv",
                "description": "Encode RNA expression - liver",
            },
            "gastrocnemius": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=gastrocnemius",
                "output": "encode.gastrocnemius.expression.tsv",
                "description": "Encode RNA expression - gastrocnemius",
            },
            "left_cerebral_cortex": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=left+cerebral+cortex",
                "output": "encode.left_cerebral_cortex.expression.tsv",
                "description": "Encode RNA expression - left cerebral cortex",
            },
            "layer_of_hippocampus": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=layer+of+hippocampus",
                "output": "encode.layer_of_hippocampus.expression.tsv",
                "description": "Encode RNA expression - layer of hippocampus",
            },
            "g1e_er4": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=G1E-ER4",
                "output": "encode.g1e_er4.expression.tsv",
                "description": "Encode RNA expression - G1E-ER4",
            },
            "forebrain": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=forebrain",
                "output": "encode.forebrain.expression.tsv",
                "description": "Encode RNA expression - forebrain",
            },
            "hindbrain": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=hindbrain",
                "output": "encode.hindbrain.expression.tsv",
                "description": "Encode RNA expression - hindbrain",
            },
            "midbrain": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=midbrain",
                "output": "encode.midbrain.expression.tsv",
                "description": "Encode RNA expression - midbrain",
            },
            "limb": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=limb",
                "output": "encode.limb.expression.tsv",
                "description": "Encode RNA expression - limb",
            },
            "kidney": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=kidney",
                "output": "encode.kidney.expression.tsv",
                "description": "Encode RNA expression - kidney",
            },
            "lung": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=lung",
                "output": "encode.lung.expression.tsv",
                "description": "Encode RNA expression - lung",
            },
            "embryonic_facial_prominence": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=embryonic+facial+prominence",
                "output": "encode.embryonic_facial_prominence.expression.tsv",
                "description": "Encode RNA expression - embryonic facial prominence",
            },
            "inflammation_experienced_regulatory_t_cells": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=inflammation-experienced+regulatory+T-cells",
                "output": "encode.inflammation_experienced_regulatory_t_cells.expression.tsv",
                "description": "Encode RNA expression - inflammation-experienced regulatory T-cells",
            },
            "neural_tube": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=neural+tube",
                "output": "encode.neural_tube.expression.tsv",
                "description": "Encode RNA expression - neural tube",
            },
            "mel": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=MEL",
                "output": "encode.mel.expression.tsv",
                "description": "Encode RNA expression - MEL",
            },
            "activated_regulatory_t_cell": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=activated+regulatory+T+cell",
                "output": "encode.activated_regulatory_t_cell.expression.tsv",
                "description": "Encode RNA expression - activated regulatory T cell",
            },
            "dendritic_cell": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=dendritic+cell",
                "output": "encode.dendritic_cell.expression.tsv",
                "description": "Encode RNA expression - dendritic cell",
            },
            "stomach": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=stomach",
                "output": "encode.stomach.expression.tsv",
                "description": "Encode RNA expression - stomach",
            },
            "g1e": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=G1E",
                "output": "encode.g1e.expression.tsv",
                "description": "Encode RNA expression - G1E",
            },
            "intestine": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=intestine",
                "output": "encode.intestine.expression.tsv",
                "description": "Encode RNA expression - intestine",
            },
            "regulatory_t_cell": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=regulatory+T+cell",
                "output": "encode.regulatory_t_cell.expression.tsv",
                "description": "Encode RNA expression - regulatory T cell",
            },
            "spleen": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=spleen",
                "output": "encode.spleen.expression.tsv",
                "description": "Encode RNA expression - spleen",
            },
            "ch12_lx": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=CH12.LX",
                "output": "encode.ch12_lx.expression.tsv",
                "description": "Encode RNA expression - CH12.LX",
            },
            "brain": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=brain",
                "output": "encode.brain.expression.tsv",
                "description": "Encode RNA expression - brain",
            },
            "central_nervous_system": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=central+nervous+system",
                "output": "encode.central_nervous_system.expression.tsv",
                "description": "Encode RNA expression - central nervous system",
            },
            "erythroblast": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=erythroblast",
                "output": "encode.erythroblast.expression.tsv",
                "description": "Encode RNA expression - erythroblast",
            },
            "megakaryocyte": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=megakaryocyte",
                "output": "encode.megakaryocyte.expression.tsv",
                "description": "Encode RNA expression - megakaryocyte",
            },
            "megakaryocyte_erythroid_progenitor_cell": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=megakaryocyte-erythroid+progenitor+cell",
                "output": "encode.megakaryocyte_erythroid_progenitor_cell.expression.tsv",
                "description": "Encode RNA expression - megakaryocyte-erythroid progenitor cell",
            },
            "small_intestine": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=small+intestine",
                "output": "encode.small_intestine.expression.tsv",
                "description": "Encode RNA expression - small intestine",
            },
            "testis": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=testis",
                "output": "encode.testis.expression.tsv",
                "description": "Encode RNA expression - testis",
            },
            "thymus": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=thymus",
                "output": "encode.thymus.expression.tsv",
                "description": "Encode RNA expression - thymus",
            },
            "c3h10t1_2": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=C3H10T1%2F2",
                "output": "encode.c3h10t1_2.expression.tsv",
                "description": "Encode RNA expression - C3H10T1/2",
            },
            "cerebellum": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=cerebellum",
                "output": "encode.cerebellum.expression.tsv",
                "description": "Encode RNA expression - cerebellum",
            },
            "common_myeloid_progenitor": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=common+myeloid+progenitor",
                "output": "encode.common_myeloid_progenitor.expression.tsv",
                "description": "Encode RNA expression - common myeloid progenitor",
            },
            "cortical_plate": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=cortical+plate",
                "output": "encode.cortical_plate.expression.tsv",
                "description": "Encode RNA expression - cortical plate",
            },
            "effector_memory_cd4_positive_alpha_beta_t_cell": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=effector+memory+CD4-positive%2C+alpha-beta+T+cell",
                "output": "encode.effector_memory_cd4_positive_alpha_beta_t_cell.expression.tsv",
                "description": "Encode RNA expression - effector memory CD4-positive, alpha-beta T cell",
            },
            "erythroid_progenitor_cell": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=erythroid+progenitor+cell",
                "output": "encode.erythroid_progenitor_cell.expression.tsv",
                "description": "Encode RNA expression - erythroid progenitor cell",
            },
            "granulocyte_monocyte_progenitor_cell": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=granulocyte+monocyte+progenitor+cell",
                "output": "encode.granulocyte_monocyte_progenitor_cell.expression.tsv",
                "description": "Encode RNA expression - granulocyte monocyte progenitor cell",
            },
            "hematopoietic_stem_cell": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=hematopoietic+stem+cell",
                "output": "encode.hematopoietic_stem_cell.expression.tsv",
                "description": "Encode RNA expression - hematopoietic stem cell",
            },
            "megakaryocyte_progenitor_cell": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=megakaryocyte+progenitor+cell",
                "output": "encode.megakaryocyte_progenitor_cell.expression.tsv",
                "description": "Encode RNA expression - megakaryocyte progenitor cell",
            },
            "ovary": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=ovary",
                "output": "encode.ovary.expression.tsv",
                "description": "Encode RNA expression - ovary",
            },
            "placenta": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=placenta",
                "output": "encode.placenta.expression.tsv",
                "description": "Encode RNA expression - placenta",
            },
            "urinary_bladder": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=urinary+bladder",
                "output": "encode.urinary_bladder.expression.tsv",
                "description": "Encode RNA expression - urinary bladder",
            },
            "c2c12": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=C2C12",
                "output": "encode.c2c12.expression.tsv",
                "description": "Encode RNA expression - C2C12",
            },
            "es_bruce4": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=ES-Bruce4",
                "output": "encode.es_bruce4.expression.tsv",
                "description": "Encode RNA expression - ES-Bruce4",
            },
            "es_e14": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=ES-E14",
                "output": "encode.es_e14.expression.tsv",
                "description": "Encode RNA expression - ES-E14",
            },
            "f121_9": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=F121-9",
                "output": "encode.f121_9.expression.tsv",
                "description": "Encode RNA expression - F121-9",
            },
            "patski": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&file.biosample_ontology.term_name=Patski",
                "output": "encode.patski.expression.tsv",
                "description": "Encode RNA expression - Patski",
            },
            "adipose_tissue": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=adipose+tissue",
                "output": "encode.adipose_tissue.expression.tsv",
                "description": "Encode RNA expression - adipose tissue",
            },
            "bone_marrow": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=bone+marrow",
                "output": "encode.bone_marrow.expression.tsv",
                "description": "Encode RNA expression - bone marrow",
            },
            "bone_marrow_macrophage": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=bone+marrow+macrophage",
                "output": "encode.bone_marrow_macrophage.expression.tsv",
                "description": "Encode RNA expression - bone marrow macrophage",
            },
            "brown_adipose_tissue": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=brown+adipose+tissue",
                "output": "encode.brown_adipose_tissue.expression.tsv",
                "description": "Encode RNA expression - brown adipose tissue",
            },
            "colon": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=colon",
                "output": "encode.colon.expression.tsv",
                "description": "Encode RNA expression - colon",
            },
            "duodenum": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=duodenum",
                "output": "encode.duodenum.expression.tsv",
                "description": "Encode RNA expression - duodenum",
            },
            "embryonic_fibroblast": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=embryonic+fibroblast",
                "output": "encode.embryonic_fibroblast.expression.tsv",
                "description": "Encode RNA expression - embryonic fibroblast",
            },
            "frontal_cortex": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=frontal+cortex",
                "output": "encode.frontal_cortex.expression.tsv",
                "description": "Encode RNA expression - frontal cortex",
            },
            "gonadal_fat_pad": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=gonadal+fat+pad",
                "output": "encode.gonadal_fat_pad.expression.tsv",
                "description": "Encode RNA expression - gonadal fat pad",
            },
            "hematopoietic_multipotent_progenitor_cell": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=hematopoietic+multipotent+progenitor+cell",
                "output": "encode.hematopoietic_multipotent_progenitor_cell.expression.tsv",
                "description": "Encode RNA expression - hematopoietic multipotent progenitor cell",
            },
            "large_intestine": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=large+intestine",
                "output": "encode.large_intestine.expression.tsv",
                "description": "Encode RNA expression - large intestine",
            },
            "leukemia_stem_cell": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=leukemia+stem+cell",
                "output": "encode.leukemia_stem_cell.expression.tsv",
                "description": "Encode RNA expression - leukemia stem cell",
            },
            "mammary_gland": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=mammary+gland",
                "output": "encode.mammary_gland.expression.tsv",
                "description": "Encode RNA expression - mammary gland",
            },
            "monocyte": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=monocyte",
                "output": "encode.monocyte.expression.tsv",
                "description": "Encode RNA expression - monocyte",
            },
            "myocyte": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=myocyte",
                "output": "encode.myocyte.expression.tsv",
                "description": "Encode RNA expression - myocyte",
            },
            "naive_thymus_derived_cd4_positive_alpha_beta_t_cell": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=naive+thymus-derived+CD4-positive%2C+alpha-beta+T+cell",
                "output": "encode.naive_thymus_derived_cd4_positive_alpha_beta_t_cell.expression.tsv",
                "description": "Encode RNA expression - naive thymus-derived CD4-positive, alpha-beta T cell",
            },
            "neutrophil": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=neutrophil",
                "output": "encode.neutrophil.expression.tsv",
                "description": "Encode RNA expression - neutrophil",
            },
            "olfactory_bulb": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=olfactory+bulb",
                "output": "encode.olfactory_bulb.expression.tsv",
                "description": "Encode RNA expression - olfactory bulb",
            },
            "pancreas": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=pancreas",
                "output": "encode.pancreas.expression.tsv",
                "description": "Encode RNA expression - pancreas",
            },
            "sigmoid_colon": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=sigmoid+colon",
                "output": "encode.sigmoid_colon.expression.tsv",
                "description": "Encode RNA expression - sigmoid colon",
            },
            "skeletal_muscle_tissue": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=skeletal+muscle+tissue",
                "output": "encode.skeletal_muscle_tissue.expression.tsv",
                "description": "Encode RNA expression - skeletal muscle tissue",
            },
            "subcutaneous_adipose_tissue": {
                "url": "https://www.encodeproject.org/report.tsv?type=RNAExpression&dataset.replicates.library.biosample.donor.organism.scientific_name=Mus+musculus&file.assay_title=total+RNA-seq&file.assay_title=polyA+plus+RNA-seq&field=%40id&field=expression.gene_id&field=expression.tpm&file.biosample_ontology.term_name=subcutaneous+adipose+tissue",
                "output": "encode.subcutaneous_adipose_tissue.expression.tsv",
                "description": "Encode RNA expression - subcutaneous adipose tissue",
            },
        },
        "fantom5": {
            "url": "https://tusco.s3.us-east-1.amazonaws.com/E-MTAB-3579-query-results.tsv",
            "output": "fantom5.tpms.tsv",
            "description": "FANTOM5 TPMs for Mouse",
        },
        "bgee": {
            "bgee_expression": {
                "url": "https://www.bgee.org/ftp/current/download/calls/expr_calls/Mus_musculus_expr_simple.tsv.gz",
                "output": "bgee.expr_simple.tsv.gz",
                "description": "Bgee Mouse Expression",
            },
            "bgee_mapping": {
                "url": "https://www.bgee.org/ftp/current/download/processed_expr_values/rna_seq/Mus_musculus/Mus_musculus_RNA-Seq_experiments_libraries.tar.gz",
                "output": "bgee.rna_seq_experiments_libraries.tar.gz",
                "description": "Bgee Mouse RNA-Seq Experiments and Libraries",
            },
        },
        "hrt_atlas": {
            "url": "https://housekeeping.unicamp.br/Housekeeping_GenesMouse.csv",
            "output": "hrt_atlas.housekeeping_genes.csv",
            "description": "HRT Atlas Housekeeping Genes for Mouse",
        },
    },
    "splicing_tss": {
        "recount3": {
            "url": "https://ftp.ebi.ac.uk/pub/databases/havana/recount3/m39_recount3_final_sorted.pass1.bed.gz",
            "output": "recount3.pass1.bed.gz",
            "description": "Recount3 Pass1 Mouse",
        },
        "refTSS": {
            "url": "https://reftss.riken.jp/datafiles/current/mouse/refTSS_v4.1_mouse_coordinate.mm39.bed.gz",
            "output": "refTSS.v4.1.mouse.mm39.bed.gz",
            "description": "RefTSS v4.1 Mouse MM39",
        },
    },
    "manual_filter": {
        "url": "https://tusco.s3.us-east-1.amazonaws.com/manual_filtered_mouse.txt",
        "output": "manual_filtered.txt",
        "description": "Manual filtering of TUSCO dataset for mouse",
    },
}
