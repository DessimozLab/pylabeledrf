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
        self.assertEqual(res,3)

        # introduce 5 random edits and compute the distance
        t2 = mutateLabeledTree(self.t1, 5)
        res = computeLRF(self.t1,t2)
        self.assertEqual(res,5)

    def test_computeLRF2(self):
        # randomise the labels and compute the distance
        t3 = randomLabels(self.t1)
        t4 = mutateLabeledTree(t3, 5)
        res = computeLRF(t3,t4)
        self.assertEqual(res,5)

if __name__ == '__main__':
        unittest.main()
