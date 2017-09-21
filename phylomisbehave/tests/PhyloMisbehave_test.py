import unittest
import os
from  phylomisbehave.PhyloMisbehave import PhyloMisbehave

class TestPhyloMisbehave(unittest.TestCase):
	
	'''Check the windows are counted correctly'''
	def test_sliding_window(self):
		
		p = PhyloMisbehave(options={'multifasta': 'abc', 'output_prefix': 'out', 'gff_file': 'gff', 'faa_file': 'faa_file', 'threads': 1, 'evalue': 0.001, 'min_snps': 5})
		p.sliding_window_size = 25
		window_count, window_start = p.sliding_window(100, [5, 10, 20, 40, 2, 99])
		self.assertEqual([4, 1, 0, 1] == window_count)
		self.assertEqual([0, 24, 49, 74] == window_start)
