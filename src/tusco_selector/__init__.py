"""
TUSCO Selector package for transcript selection and processing.
"""

__version__ = "0.1.0"

from .cli import list_species, run
from .pipeline.introverse_filter import (
    read_csv_and_get_ensembl_ids as read_introverse_csv,
)
