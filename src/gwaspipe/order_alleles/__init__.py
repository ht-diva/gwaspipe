"""
Order alleles module for GWASPipe.

This module provides functionality for ordering alleles in summary statistics
data according to custom sorting rules and updating statistics accordingly.
"""

from .constants import ORDER_MAPPING, TRANSLATE_TABLE_ORDER
from .main import order_alleles
from .snpid import build_snpids, parallelbuildsnpid
from .sorting import custom_alleles_sort
from .vectorized import (
    _orderalleles_status_vec,
    orderalleles_status,
    parallelorderalleles_status,
    vectorizedorderalleles_status,
)

__all__ = [
    "custom_alleles_sort",
    "vectorizedorderalleles_status",
    "orderalleles_status",
    "parallelorderalleles_status",
    "build_snpids",
    "parallelbuildsnpid",
    "order_alleles",
    "_orderalleles_status_vec",
    "ORDER_MAPPING",
    "TRANSLATE_TABLE_ORDER",
]
