"""Ensembl REST API service integration."""

from __future__ import annotations

import json
import logging
import time
from typing import List, Set
from urllib.parse import quote

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from tusco_selector.logging_utils import get_logger

logger = get_logger(__name__)

# Ensembl REST API configuration
ENSEMBL_REST_URL = "https://rest.ensembl.org"
ENSEMBL_BATCH_SIZE = 200  # Max IDs per batch request
ENSEMBL_RATE_LIMIT_DELAY = 0.1  # Seconds between requests


def _create_session() -> requests.Session:
    """Create a requests session with retry logic.
    
    Returns:
        Configured requests Session
    """
    session = requests.Session()
    
    # Configure retry strategy
    retry_strategy = Retry(
        total=3,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"],
    )
    
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    
    # Set headers
    session.headers.update({
        "Content-Type": "application/json",
        "Accept": "application/json",
    })
    
    return session


def lookup_region_genes(location: str, species: str = "human") -> Set[str]:
    """Look up genes in a genomic region using Ensembl REST API.
    
    Args:
        location: Genomic location in format "chr:start-end" (e.g., "1:1000000-2000000")
        species: Species name for Ensembl (default: "human")
        
    Returns:
        Set of gene IDs found in the region
    """
    logger.info(f"Looking up genes in region: {location}")
    
    # Parse location
    try:
        chrom, coords = location.split(":")
        start, end = coords.split("-")
        start = int(start)
        end = int(end)
    except (ValueError, IndexError) as e:
        logger.error(f"Invalid location format '{location}': {e}")
        return set()
    
    # Build API URL
    url = f"{ENSEMBL_REST_URL}/overlap/region/{species}/{chrom}:{start}-{end}"
    params = {
        "feature": "gene",
        "content-type": "application/json",
    }
    
    session = _create_session()
    
    try:
        # Make request
        response = session.get(url, params=params, timeout=30)
        response.raise_for_status()
        
        # Parse response
        data = response.json()
        
        # Extract gene IDs
        gene_ids = set()
        for feature in data:
            if feature.get("biotype") == "protein_coding":
                # Try different ID fields
                gene_id = (
                    feature.get("gene_id") or
                    feature.get("id") or
                    feature.get("external_name")
                )
                if gene_id:
                    gene_ids.add(gene_id)
        
        logger.info(f"Found {len(gene_ids)} protein-coding genes in region {location}")
        return gene_ids
        
    except requests.exceptions.RequestException as e:
        logger.error(f"Error querying Ensembl for region {location}: {e}")
        return set()
    except (json.JSONDecodeError, KeyError) as e:
        logger.error(f"Error parsing Ensembl response: {e}")
        return set()
    finally:
        session.close()
        # Rate limiting
        time.sleep(ENSEMBL_RATE_LIMIT_DELAY)


def transcripts_to_genes(
    transcript_ids: List[str],
    species: str = "human",
    batch_size: int = ENSEMBL_BATCH_SIZE,
) -> Set[str]:
    """Convert transcript IDs to gene IDs using Ensembl REST API.
    
    Handles batching and retry logic for large lists of transcripts.
    
    Args:
        transcript_ids: List of transcript IDs to convert
        species: Species name for Ensembl (default: "human")
        batch_size: Number of IDs to process per request
        
    Returns:
        Set of corresponding gene IDs
    """
    if not transcript_ids:
        return set()
    
    logger.info(f"Converting {len(transcript_ids)} transcript IDs to gene IDs")
    
    # Remove duplicates
    unique_transcripts = list(set(transcript_ids))
    
    gene_ids = set()
    session = _create_session()
    
    try:
        # Process in batches
        for i in range(0, len(unique_transcripts), batch_size):
            batch = unique_transcripts[i:i + batch_size]
            logger.debug(f"Processing batch {i//batch_size + 1} with {len(batch)} transcripts")
            
            # Build POST request
            url = f"{ENSEMBL_REST_URL}/lookup/id/{species}"
            data = {"ids": batch}
            
            try:
                # Make request
                response = session.post(
                    url,
                    json=data,
                    timeout=60,
                )
                response.raise_for_status()
                
                # Parse response
                result = response.json()
                
                # Extract gene IDs
                for transcript_id, info in result.items():
                    if info and isinstance(info, dict):
                        gene_id = info.get("Parent") or info.get("gene_id")
                        if gene_id:
                            gene_ids.add(gene_id)
                
            except requests.exceptions.RequestException as e:
                logger.error(f"Error in batch {i//batch_size + 1}: {e}")
                # Continue with next batch rather than failing completely
                continue
            except (json.JSONDecodeError, KeyError) as e:
                logger.error(f"Error parsing response for batch {i//batch_size + 1}: {e}")
                continue
            finally:
                # Rate limiting between batches
                time.sleep(ENSEMBL_RATE_LIMIT_DELAY)
        
        logger.info(f"Converted to {len(gene_ids)} unique gene IDs")
        return gene_ids
        
    finally:
        session.close()


def get_gene_info(gene_id: str, species: str = "human") -> dict:
    """Get detailed information about a gene from Ensembl.
    
    Args:
        gene_id: Ensembl gene ID
        species: Species name for Ensembl
        
    Returns:
        Dictionary with gene information
    """
    logger.debug(f"Getting info for gene: {gene_id}")
    
    url = f"{ENSEMBL_REST_URL}/lookup/id/{gene_id}"
    params = {"expand": 1}  # Include additional info
    
    session = _create_session()
    
    try:
        response = session.get(url, params=params, timeout=30)
        response.raise_for_status()
        return response.json()
        
    except requests.exceptions.RequestException as e:
        logger.error(f"Error getting info for gene {gene_id}: {e}")
        return {}
    finally:
        session.close()
        time.sleep(ENSEMBL_RATE_LIMIT_DELAY)


__all__ = [
    'lookup_region_genes',
    'transcripts_to_genes',
    'get_gene_info',
]