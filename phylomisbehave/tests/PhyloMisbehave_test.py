import unittest
import os
from  phylomisbehave.PhyloMisbehave import PhyloMisbehave

class Options:
	def __init__(self, multifasta, output_prefix, gff_file, faa_file, threads, evalue, min_snps,verbose):	
		self.multifasta    = multifasta 
		self.output_prefix = output_prefix
		self.gff_file      = gff_file
		self.faa_file      = faa_file
		self.threads       = threads
		self.evalue        = evalue
		self.min_snps      = min_snps
		self.verbose       = verbose

class TestPhyloMisbehave(unittest.TestCase):
	
	'''Check the windows are counted correctly'''
	def test_sliding_window(self):
		
		options = Options('mfa','out','gff','faa',1,0.001,5, False)
		p = PhyloMisbehave(options)
		p.sliding_window_size = 25
		window_count, window_start = p.sliding_window(100, [5, 10, 20, 40, 2, 78])
		self.assertEqual([4, 1, 0, 1],  window_count)
		self.assertEqual([0, 24, 49, 74], window_start)
