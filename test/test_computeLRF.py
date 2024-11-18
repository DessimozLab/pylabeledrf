import unittest
from pylabeledrf.computeLRF import *
import dendropy
import sys
import os

TESTDATA = os.path.join(os.path.dirname(__file__), 'p53.nhx')

class UtilTest(unittest.TestCase):
    def setUp(self):
        self.taxa = dendropy.TaxonNamespace()
        self.p53 = dendropy.Tree.get_from_path(TESTDATA, 
            'newick', taxon_namespace=self.taxa)
        self.t1 = parseEnsemblLabels(self.p53)

    def test_computeLRF1(self):
        # introduce 3 random edits and compute the distance
        t2 = mutateLabeledTree(self.t1, 3)
        res = computeLRF(self.t1,t2)
        self.assertLessEqual(res,3)

        # introduce 5 random edits and compute the distance
        t2 = mutateLabeledTree(self.t1, 5)
        res = computeLRF(self.t1,t2)
        self.assertLessEqual(res,5)

    def test_computeLRF2(self):
        # randomise the labels and compute the distance
        t3 = randomLabels(self.t1)
        t4 = mutateLabeledTree(t3, 5)
        res = computeLRF(t3,t4)
        self.assertLessEqual(res,5)


class ErrorTests(unittest.TestCase):
    def test_raises_error_if_not_same_taxon_namespace(self):
        t1 = dendropy.Tree.get_from_string("(A,(B,C)duplication);", schema="newick")
        t2 = dendropy.Tree.get_from_string("(A,(B,C)duplication);", schema="newick")
        with self.assertRaises(ValueError):
            computeLRF(t1, t2)

    def test_raises_no_input(self):
        with self.assertRaises(TypeError):
            computeLRF()

    def test_raises_invalid_argument_types(self):
        with self.assertRaises(ValueError):
            computeLRF('((A,B),C);', '(A,(B,C));')


if __name__ == '__main__':
        unittest.main()
