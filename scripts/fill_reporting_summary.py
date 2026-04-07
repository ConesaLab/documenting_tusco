#!/usr/bin/env python3
"""Fill the Nature Portfolio Reporting Summary XFA form for the TUSCO paper.

Reads nr-reporting-summary.pdf, populates fields via XFA datasets XML,
and writes nr-reporting-summary-filled.pdf.

Usage:
    python scripts/fill_reporting_summary.py
"""

import xml.etree.ElementTree as ET
import zlib
from pathlib import Path

import pypdf
from pypdf.generic import NameObject, NumberObject

REPO = Path(__file__).resolve().parent.parent
INPUT_PDF = REPO / "submission" / "templates" / "nr-reporting-summary.pdf"
OUTPUT_PDF = REPO / "submission" / "01_for_upload" / "Reporting_Summary.pdf"

# ── XFA namespace ──────────────────────────────────────────────────────
NS = {"xfa": "http://www.xfa.org/schema/xfa-data/1.0/"}


# ── Answers ────────────────────────────────────────────────────────────
# Radio/checkbox values: 1 = N/A, 2 = Yes/confirmed

HEADER = {
    "author": "Ana Conesa",
    "lastupdate": "2026-03-19",
}

# SelectField > Type: 1=Life sciences, 2=Behavioural, 3=EEE
FIELD_SELECTION = "1"

# Statistics confirmation checkboxes (under Statistics > questions)
STATISTICS_QUESTIONS = {
    "samplesize": "2",    # confirmed
    "collection": "2",    # confirmed
    "statstests": "2",    # confirmed
    "covariates": "2",    # confirmed
    "assumptions": "2",   # confirmed
    "descstats": "2",     # confirmed
    "nullhypothesis": "1",  # N/A — permutation test has no classical test statistic/df/CI
    "bayesian": "1",        # N/A
    "hierarchical": "1",    # N/A
    "effectsizes": "1",     # N/A
}

SOFTWARE = {
    "collectioninfo": (
        "Long-read RNA-seq data were generated in-house using PacBio "
        "Sequel IIe (mouse kidney/brain, Iso-Seq SMRTbell prep kit 3.0) "
        "and obtained from public repositories: LRGASP Consortium (ONT "
        "MinION R10.4.1, GridION R9.4.1, PacBio Sequel II) via Synapse "
        "(syn25007472), RNA degradation ONT cDNA data via ENA "
        "(PRJEB53210) and EGA (EGAS00001006542). Short-read Illumina "
        "data were aligned with STAR for junction support."
    ),
    "analysisinfo": (
        "R (v4.3.3) with ggplot2 (v3.5.1), Python 3.11, "
        "SQANTI3 (v5.4), minimap2 (v2.28), FLAIR (v2.2.0), "
        "STAR (v2.7.11b), Bambu (v3.8.3), StringTie2 (v2.2.3), "
        "PacBio Iso-Seq pipeline (v4.3.0), "
        "MAJIQ-L (v3.0.16), AlphaGenome (v0.1.0), "
        "samtools (v1.21), Ensembl BioMart. "
        "All analysis scripts available at "
        "https://github.com/ConesaLab/documenting_tusco"
    ),
}

DATA_AVAILABILITY = (
    "Publicly available LRS datasets used in this study include the "
    "LRGASP Consortium data for human WTC11 iPSC and mouse ES cells, "
    "accessible at the Synapse repository under accession "
    "syn25007472/wiki/608702. The RNA degradation datasets are deposited "
    "in the European Nucleotide Archive (ENA) under accession PRJEB53210, "
    "with supplementary FASTQ data for post-mortem brain samples available "
    "from the European Genome-Phenome Archive (EGA) under accession "
    "EGAS00001006542. The Replicated Mouse Brain and Kidney dataset, "
    "generated in-house, has been submitted in its entirety to the ENA "
    'under accessions PRJEB85167 (secondary accession ERP168594, "Sequencing '
    'of mice brains using Iso-Seq") and PRJEB94912 (secondary accession '
    'ERP177675, "Sequencing of mice kidneys using Iso-Seq").\n'
    "\n"
    "For benchmarking, the universal and tissue-specific TUSCO gene sets "
    "for human and mouse can be freely obtained at "
    "https://github.com/ConesaLab/SQANTI3 and at the dedicated TUSCO "
    "benchmarking framework portal https://tusco.uv.es. An example of "
    "the TUSCO report is available at: "
    "https://github.com/ConesaLab/SQANTI3/blob/master/example/"
    "tusco_output/tusco_report.html. Source data are provided with this "
    "paper. All data are also available from the corresponding author "
    "upon request."
)

# Life Sciences Study Design
LIFE_SCIENCES = {
    "howselectedsamplesize": (
        "No sample-size calculation was performed. TUSCO gene set sizes "
        "(human universal: 46 genes, mouse universal: 32 genes) resulted "
        "from applying pre-established multi-database filtering criteria. "
        "Mouse kidney/brain: 5 biological replicates per tissue, chosen to "
        "enable replicate-combination consensus analysis. LRGASP: 6 pipeline "
        "combinations per species (defined by consortium). RNA degradation: "
        "17 ONT cDNA samples (RIN 7.2\u20139.9)."
    ),
    "exclusions": (
        "Genes were excluded from the TUSCO gene set based on "
        "pre-established criteria: failure to achieve cross-annotation "
        "single-isoform agreement (GENCODE, RefSeq, MANE/Ensembl), "
        "insufficient expression prevalence (<95% of tissues for human, "
        "<90% for mouse), splice invariance filter failure, or TSS "
        "inconsistency. No data points were excluded from downstream "
        "analyses after gene set selection."
    ),
    "limitations": (
        "Mouse kidney/brain: 5 biological replicates per tissue; all "
        "possible replicate combinations computed with mean and 95% CI "
        "reported. Cross-platform (PacBio/ONT) and cross-library "
        "(dRNA/cDNA/hybrid) validation performed. TUSCO benchmarks "
        "validated against SIRVs (cosine similarity 0.9493\u20130.9992) and "
        "MAJIQ-L junction-level analysis. All findings were consistent "
        "across species, platforms, and pipelines."
    ),
    "methodrandomization": (
        "Not applicable. This is a computational/bioinformatic benchmarking "
        "study with no experimental group allocation requiring "
        "randomization."
    ),
    "blindinggroupallocation": (
        "Not applicable. Computational analyses were performed using "
        "pre-defined pipelines and criteria; no subjective assessment "
        "requiring blinding was involved."
    ),
}

# Materials & Experimental Systems checklist: 1=N/A, 2=Yes
MATERIALS_CHECKLIST = {
    "antibodies": "1",
    "celllines": "1",
    "palaeontology": "1",
    "animals": "2",       # Yes — C57BL/6J mice
    "clinical": "1",
    "durc": "1",
    "plants": "1",
}

# Methods checklist: 1=N/A, 2=Yes
METHODS_CHECKLIST = {
    "chip-seq": "1",
    "flowcytometry": "1",
    "mri": "1",
}

# Animals and Other Research Organisms
ANIMALS = {
    "animaldetails": (
        "C57BL/6J wild-type mice, male, 3 months old. Five animals used "
        "for kidney tissue and the same cohort for brain tissue. Tissues "
        "homogenized using FastPrep; RNA extracted with Maxwell 16 LEV "
        "simplyRNA Purification Kit."
    ),
    "wildanimals": "Not applicable.",
    # First <ethics> in Animals = Reporting on sex
    "ethics_sex": (
        "Only male mice were used. The study focused on transcriptome "
        "benchmarking methodology rather than sex-specific biology; a "
        "single sex was used to reduce biological variability."
    ),
    "fieldsamples": "Not applicable.",
    # Second <ethics> in Animals = Ethics oversight
    "ethics_oversight": (
        "All animal handling and experimental procedures were conducted "
        "in accordance with European Union 86/609/EEC and Spanish "
        "RD1201/2005 guidelines, following protocols approved by the "
        "Ethics Committee on Experimental Research of the University "
        "of Valencia."
    ),
}


# ── Helper: set text of an XML element, creating it if needed ──────────
def set_text(parent, tag, value):
    """Set the text content of a child element, creating it if absent."""
    el = parent.find(tag)
    if el is None:
        el = ET.SubElement(parent, tag)
    el.text = value


def find_or_create(parent, tag):
    """Find a child element or create it."""
    el = parent.find(tag)
    if el is None:
        el = ET.SubElement(parent, tag)
    return el


# ── Main ───────────────────────────────────────────────────────────────
def fill_form():
    reader = pypdf.PdfReader(str(INPUT_PDF))
    writer = pypdf.PdfWriter(clone_from=reader)

    # Get XFA array from writer's AcroForm
    acroform = writer._root_object["/AcroForm"]
    xfa_array = acroform["/XFA"]

    # The datasets stream is at index 9 in the XFA array
    datasets_obj = xfa_array[9].get_object()
    datasets_xml = datasets_obj.get_data().decode("utf-8")

    # Parse the datasets XML
    root = ET.fromstring(datasets_xml)

    # Navigate: root = <xfa:datasets>, child = <xfa:data>, grandchild = <form>
    xfa_data = root.find("xfa:data", NS)
    form = xfa_data.find("form")
    main = form.find("Maincontent")

    # ── Header ─────────────────────────────────────────────────────
    header = find_or_create(main, "Header")
    set_text(header, "author", HEADER["author"])
    set_text(header, "lastupdate", HEADER["lastupdate"])

    # ── Statistics questions ───────────────────────────────────────
    universal = find_or_create(main, "Universal")
    stats = find_or_create(universal, "Statistics")
    questions = find_or_create(stats, "questions")
    for qname, qval in STATISTICS_QUESTIONS.items():
        set_text(questions, qname, qval)

    # ── Software ───────────────────────────────────────────────────
    software = find_or_create(universal, "Software")
    set_text(software, "collectioninfo", SOFTWARE["collectioninfo"])
    set_text(software, "analysisinfo", SOFTWARE["analysisinfo"])

    # ── Data Availability ──────────────────────────────────────────
    dataavail = find_or_create(universal, "Dataavail")
    das = find_or_create(dataavail, "DAS")
    set_text(das, "howselectedsamplesize", DATA_AVAILABILITY)

    # ── Humans (not applicable) ───────────────────────────────────
    humans = find_or_create(main, "Humans")
    human_details = humans.findall("humandetails")
    for hd in human_details:
        if not hd.text:
            hd.text = "NA"
    set_text(humans, "recruitment", "NA")
    # Humans > ethics = ethics oversight for human participants
    humans_ethics = humans.find("ethics")
    if humans_ethics is not None and not humans_ethics.text:
        humans_ethics.text = "NA"

    # ── SelectField (Life sciences = 1) ────────────────────────────
    selectfield = find_or_create(main, "SelectField")
    set_text(selectfield, "Type", FIELD_SELECTION)

    # ── Life Sciences Study Design ─────────────────────────────────
    life = find_or_create(main, "Life")
    for fname, fval in LIFE_SCIENCES.items():
        set_text(life, fname, fval)

    # ── Materials & Experimental Systems checklist ─────────────────
    method_modules = find_or_create(main, "MethodModules")
    msr = find_or_create(method_modules, "MethodspecificReporting")
    materials = find_or_create(msr, "materials")
    for mname, mval in MATERIALS_CHECKLIST.items():
        set_text(materials, mname, mval)

    # ── Methods checklist ──────────────────────────────────────────
    methods = find_or_create(msr, "methods")
    for mname, mval in METHODS_CHECKLIST.items():
        set_text(methods, mname, mval)

    # ── Animals section ────────────────────────────────────────────
    animals = find_or_create(main, "Animals")
    set_text(animals, "animaldetails", ANIMALS["animaldetails"])
    set_text(animals, "wildanimals", ANIMALS["wildanimals"])
    set_text(animals, "fieldsamples", ANIMALS["fieldsamples"])

    # Handle the two <ethics> elements in Animals:
    # First = Reporting on sex, Second = Ethics oversight
    ethics_elements = animals.findall("ethics")
    if len(ethics_elements) >= 2:
        ethics_elements[0].text = ANIMALS["ethics_sex"]
        ethics_elements[1].text = ANIMALS["ethics_oversight"]
    elif len(ethics_elements) == 1:
        ethics_elements[0].text = ANIMALS["ethics_sex"]
        ethics2 = ET.SubElement(animals, "ethics")
        ethics2.text = ANIMALS["ethics_oversight"]
    else:
        ethics1 = ET.SubElement(animals, "ethics")
        ethics1.text = ANIMALS["ethics_sex"]
        ethics2 = ET.SubElement(animals, "ethics")
        ethics2.text = ANIMALS["ethics_oversight"]

    # ── Serialize modified XML back ────────────────────────────────
    new_xml = ET.tostring(root, encoding="unicode", xml_declaration=False)
    # Restore the XFA namespace declaration that ET may mangle
    new_xml = new_xml.replace(
        "ns0:datasets", "xfa:datasets"
    ).replace(
        "xmlns:ns0=", "xmlns:xfa="
    ).replace(
        "ns0:data", "xfa:data"
    ).replace(
        'ns0:dataNode="dataGroup"', 'xfa:dataNode="dataGroup"'
    )

    # ── Write modified datasets back into the PDF ─────────────────
    # Compress with FlateDecode to match the original stream format
    raw_bytes = new_xml.encode("utf-8")
    compressed = zlib.compress(raw_bytes)

    datasets_ref = xfa_array[9]
    obj = (
        writer.get_object(datasets_ref)
        if isinstance(datasets_ref, pypdf.generic.IndirectObject)
        else datasets_obj
    )
    # Replace stream content: write compressed data with FlateDecode filter
    obj._data = compressed
    obj[NameObject("/Filter")] = NameObject("/FlateDecode")
    obj[NameObject("/Length")] = NumberObject(len(compressed))

    writer.write(str(OUTPUT_PDF))
    print(f"Wrote filled form to {OUTPUT_PDF}")


if __name__ == "__main__":
    fill_form()
