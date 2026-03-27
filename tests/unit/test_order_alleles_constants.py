import unittest

import pandas as pd
from gwaslab.info.g_Log import Log

from gwaspipe.order_alleles.constants import ORDER_MAPPING, TRANSLATE_TABLE_ORDER


class TestOrderMapping(unittest.TestCase):
    """Tests for ORDER_MAPPING constant."""

    def test_order_mapping_keys(self):
        """Test ORDER_MAPPING contains the four nucleotide bases."""
        self.assertEqual(set(ORDER_MAPPING.keys()), {"A", "C", "G", "T"})

    def test_order_mapping_values_are_not_null(self):
        """Test ORDER_MAPPING values are not chr(0)."""
        for v in ORDER_MAPPING.values():
            self.assertNotEqual(v, chr(0))

    def test_order_mapping_values_are_unique(self):
        """Test ORDER_MAPPING values are unique."""
        values = list(ORDER_MAPPING.values())
        self.assertEqual(len(values), len(set(values)))

    def test_order_mapping_is_sorted(self):
        """Test ORDER_MAPPING is sorted alphabetically."""
        sorted_keys = sorted(ORDER_MAPPING.keys())
        self.assertEqual(list(ORDER_MAPPING.keys()), sorted_keys)


class TestTranslateTableOrder(unittest.TestCase):
    """Tests for TRANSLATE_TABLE_ORDER constant."""

    def test_translate_table_order_exists(self):
        """Test TRANSLATE_TABLE_ORDER is created."""
        self.assertIsNotNone(TRANSLATE_TABLE_ORDER)

    def test_translate_table_order_is_bytes(self):
        """Test TRANSLATE_TABLE_ORDER is a bytes object."""
        self.assertIsInstance(TRANSLATE_TABLE_ORDER, bytes)

    def test_translate_table_order_has_correct_length(self):
        """Test TRANSLATE_TABLE_ORDER has correct length (256)."""
        self.assertEqual(len(TRANSLATE_TABLE_ORDER), 256)


if __name__ == "__main__":
    unittest.main()
