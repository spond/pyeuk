import os
import sys
import unittest
import warnings

class TestDeprecationShim(unittest.TestCase):

    def test_cyclospora_pyeuk_deprecation_warning(self):
        # Remove cached imports if any
        for key in list(sys.modules.keys()):
            if key == "cyclospora_pyeuk" or key.startswith("cyclospora_pyeuk."):
                del sys.modules[key]

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            import cyclospora_pyeuk
            self.assertTrue(any(issubclass(item.category, DeprecationWarning) for item in w))
            self.assertTrue(any("deprecated" in str(item.message) for item in w))

    def test_cyclospora_pyeuk_symbol_parity(self):
        import pyeuk
        import cyclospora_pyeuk

        self.assertEqual(cyclospora_pyeuk.__version__, pyeuk.__version__)
        self.assertIs(cyclospora_pyeuk.PyEukDistanceEngine, pyeuk.PyEukDistanceEngine)
        self.assertIs(cyclospora_pyeuk.CyclosporaClusterFinder, pyeuk.CyclosporaClusterFinder)
        self.assertIs(cyclospora_pyeuk.generate_haplotype_sheet, pyeuk.generate_haplotype_sheet)
        self.assertIs(cyclospora_pyeuk.name_haplotype, pyeuk.name_haplotype)
        self.assertIs(cyclospora_pyeuk.parse_locus_name, pyeuk.parse_locus_name)

    def test_submodule_shims(self):
        from cyclospora_pyeuk.cli import main as old_main
        from pyeuk.cli import main as new_main
        self.assertIs(old_main, new_main)

        from cyclospora_pyeuk.distance_engine import PyEukDistanceEngine as OldEngine
        from pyeuk.distance_engine import PyEukDistanceEngine as NewEngine
        self.assertIs(OldEngine, NewEngine)

        from cyclospora_pyeuk.amplicon import build_sheet as old_build_sheet
        from pyeuk.amplicon import build_sheet as new_build_sheet
        self.assertIs(old_build_sheet, new_build_sheet)


if __name__ == "__main__":
    unittest.main()
