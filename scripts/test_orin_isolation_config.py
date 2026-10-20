#!/usr/bin/env python3
import unittest
from configure_orin_localization_isolation import proposed


class IsolationConfigTest(unittest.TestCase):
    def test_primary_only_and_idempotent(self):
        original = "DEFAULT primary\nLABEL primary\n APPEND ${cbootargs} root=/dev/nvme0n1p1\nLABEL backup\n APPEND old\n"
        result = proposed(original)
        self.assertIn("${cbootargs} root=/dev/nvme0n1p1", result)
        self.assertIn("isolcpus=domain,managed_irq,10-11 irqaffinity=0-9", result)
        self.assertTrue(result.endswith("LABEL backup\n APPEND old\n"))
        self.assertEqual(proposed(result), result)

    def test_reject_conflicting_or_ambiguous_boot_config(self):
        for original in (
            "LABEL primary\n APPEND root=x isolcpus=2-3\n",
            "LABEL primary\n APPEND root=x irqaffinity=0-11\n",
            "LABEL other\n APPEND root=x\n",
            "LABEL primary\n APPEND root=x\n APPEND root=y\n",
        ):
            with self.assertRaises(ValueError):
                proposed(original)


if __name__ == "__main__":
    unittest.main()
