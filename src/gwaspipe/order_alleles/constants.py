"""
Constants for order_alleles module.
"""

from gwaslab.bd_common_data import _maketrans

# Create translation table for fast vectorized ordering
ORDER_MAPPING = {el: chr(i + 1) for i, el in enumerate(sorted(["A", "T", "C", "G"], reverse=False))}
assert all(value != chr(0) for value in ORDER_MAPPING.values()), "Mapping should not equal chr(0)"
TRANSLATE_TABLE_ORDER = _maketrans(ORDER_MAPPING)
