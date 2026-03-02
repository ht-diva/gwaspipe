"""
Order alleles module for GWASPipe.

This module provides functionality for ordering alleles in summary statistics
data according to custom sorting rules and updating statistics accordingly.
"""

from .sorting import custom_alleles_sort
from .vectorized import (
    vectorizedorderalleles_status,
    orderalleles_status,
    parallelorderalleles_status,
    _orderalleles_status_vec,
)
from .snpid import build_snpids, parallelbuildsnpid

# from .stats import _flip_allele_statistics
from .constants import ORDER_MAPPING, TRANSLATE_TABLE_ORDER
from .main import order_alleles

__all__ = [
    "custom_alleles_sort",
    "vectorizedorderalleles_status",
    "orderalleles_status",
    "parallelorderalleles_status",
    "build_snpids",
    "parallelbuildsnpid",
    "order_alleles",
    "_orderalleles_status_vec",
    # "_flip_allele_statistics",
    "ORDER_MAPPING",
    "TRANSLATE_TABLE_ORDER",
]
