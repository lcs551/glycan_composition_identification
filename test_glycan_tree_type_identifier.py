import unittest
from glycan_tree_type_identifier import (get_unique_sugars, get_sugar_order, get_linkages, organise_linkages,
                                         branches_to_sugars, check_type)

class GlycanTreeTypeIdentifierTest(unittest.TestCase):
    def setUp(self):
        self.wurcs = "WURCS=2.0/3,7,6/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]/1-1-2-3-3-3-3/a4-b1_b4-c1_c3-d1_c6-e1_e3-f1_f2-g1"

    def test_get_unique_sugars(self):
        result = get_unique_sugars(self.wurcs)
        self.assertListEqual(result, ['NAG', 'BMA', 'MAN'])

    def test_get_sugar_order(self):
        result = get_sugar_order(self.wurcs)
        self.assertListEqual(result, ['1', '1', '2', '3', '3', '3', '3'])

    def test_get_linkages(self):
        result = get_linkages(self.wurcs)
        self.assertListEqual(result, ['a4-b1', 'b4-c1', 'c3-d1', 'c6-e1', 'e3-f1', 'f2-g1'])

    def test_organise_linkages(self):
        linkages = ['a4-b1', 'b4-c1', 'c3-d1', 'c6-e1', 'e3-f1', 'f2-g1']
        result = organise_linkages(linkages)
        self.assertListEqual(result, [['d'], ['e', 'f', 'g']])

    def test_branches_to_sugars(self):
        branches = [['A', 'B', 'C']]
        sugar_alphabet_map = {'A': 'SugarA', 'B': 'SugarB', 'C': 'SugarC'}
        result = branches_to_sugars(branches, sugar_alphabet_map)
        self.assertListEqual(result, [['SugarA', 'SugarB', 'SugarC']])

    def test_check_type(self):
        self.assertEqual(check_type(self.wurcs), 'High Mannose')


if __name__ == "__main__":
    unittest.main()