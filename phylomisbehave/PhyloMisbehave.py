import matplotlib as mpl
mpl.use('Agg') # stops the need for X11
import matplotlib.pyplot as plt

import pyximport
from ete3 import Tree
import os
import subprocess
import logging
import pkg_resources
import seaborn as sns
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SearchIO
pyximport.install()

from phylomisbehave.Homoplasies import find_homoplasious_sites
import findSNPs



#todo: list which accessions are affected by clusters/homoplasy, check if they are found in monopyletic group
#todo: add argparse
#todo: add option to use vcf as input


class PhyloMisbehave:
	def __init__(self,options):
		self.multifasta = options.multifasta
		self.output_prefix = options.output_prefix
		self.gff_file = options.gff_file
		self.faa_file = options.faa_file
		
		# Todo: make this kind of thing an option
		self.sliding_window_size = 1000
		
		self.verbose = options.verbose
		self.logger = logging.getLogger(__name__)
		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)
			
		self.fasttree_exec = self.choose_executable(['FastTree','fasttree'])
		
	def choose_executable(self,list_of_executables):
		for executable in list_of_executables:
			if self.which(executable) != None:
				return executable
			  
		return ""

	def which(self, program):
		executable = program.split(" ")
		program = executable[0]
		def is_exe(fpath):
			return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
		fpath, fname = os.path.split(program)
		if fpath:
			if is_exe(program):
				return program
		else:
			for path in os.environ["PATH"].split(os.pathsep):
				exe_file = os.path.join(path, program)
				if is_exe(exe_file):
					return exe_file
				
		return None

	def sliding_window(self, genome_size, positions):
		window_count = []
		window_start = []
		#self.logger.warning('starting sliding window')
		for i in range(0, genome_size - (self.sliding_window_size -1), self.sliding_window_size):
			window_start.append(i)
			window_count.append(len(set(range(i, i+(self.sliding_window_size -1)))& set(positions)))
		return window_count, window_start
	
	
	def findFeatures(self, gff, positions):
		geneswithSNPs=[]
		with open(gff, 'r') as featurefile:
			for line in featurefile:
				if (len(line.split('\t')) > 5) and line.split('\t')[2] == 'CDS':
					SNPS = len(set(positions) & set(range(int(line.split('\t')[3])-1, int(line.split('\t')[4]))))
					if SNPS > 0:
						geneswithSNPs.append(line.split('\t')[8].split(';')[0].lstrip('ID='))
		return geneswithSNPs
	
	def featuresdict(self, gff):
		GeneDict={}
		with open(gff, 'r') as featurefile:
			for line in featurefile:
				if (len(line.split('\t')) > 5) and line.split('\t')[2] == 'CDS':
					GeneDict[line.split('\t')[8].split(';')[0].lstrip('ID=')] = (int(line.split('\t')[3])-1,\
							int(line.split('\t')[4])-1)
		return GeneDict
	
	
	def run(self):
		alignment = findSNPs.parse_fasta(open(self.multifasta))
		ref = next(alignment)[1]
		n = len(ref)
		alignment = findSNPs.parse_fasta(open(self.multifasta))
		
		gff = self.gff_file
		varsites = findSNPs.findSNPs(ref, alignment)
		alignment = findSNPs.parse_fasta(open(self.multifasta))
		
		snps_ordered, id_snps = findSNPs.setupSNPdict(varsites, alignment)
		positions = [i for i in snps_ordered]
		window_count, window_start = self.sliding_window(n, positions)
		highsnps = [(window_start[i],j) for i,j in enumerate(window_count) if j > 5]
		filtered_positions = set(positions)
		
		for i in highsnps:
			filtered_positions = filtered_positions - set(range(i[0], i[0]+1000))
		alignment = findSNPs.parse_fasta(open(self.multifasta))
		snps_ordered_filt, id_snps_filt = findSNPs.setupSNPdict(sorted(filtered_positions), alignment)
		with open(str(self.output_prefix)+'.unfiltered.fasta', 'w') as fp:
			for i in id_snps:
				if id_snps[i].count('-') > round(len(id_snps[i])*0.45):
					self.logger.warning(str(i)+' not in peak fasta!')
				else:
					fp.write('>'+ i +'\n')
					fp.write(''.join(id_snps[i])+'\n')
		subprocess.call([self.fasttree_exec+' -nt '+str(self.output_prefix)+'.unfiltered.fasta > '+str(self.output_prefix)+'.unfiltered.nwk'], stdout=subprocess.DEVNULL, shell=True)
		with open(str(self.output_prefix)+'.peakfiltered.fasta', 'w') as fp:
			for i in id_snps_filt:
				if id_snps_filt[i].count('-') > round(len(id_snps_filt[i])*0.45):
					self.logger.warning(str(i)+' not in fasta!')
				else:
					fp.write('>'+ i +'\n')
					fp.write(''.join(id_snps_filt[i])+'\n')
		subprocess.call([self.fasttree_exec+' -nt '+str(self.output_prefix)+'.peakfiltered.fasta > '+str(self.output_prefix)+'.peakfiltered.nwk'], stdout=subprocess.DEVNULL, shell=True)
		treeETE = Tree(str(self.output_prefix)+'.unfiltered.nwk')
		alignment = [(x,y) for (x,y) in findSNPs.parse_fasta(open(str(self.output_prefix)+'.unfiltered.fasta'))]
		nonhomoplasy, homoplasy = find_homoplasious_sites(alignment, treeETE)
		homoplasious_sites = [positions[x] for x in homoplasy]
		window_hom, window_start = self.sliding_window(n, homoplasious_sites)
		nonhomoplasious_sites = list(set(positions)-set(homoplasious_sites))
		alignment = findSNPs.parse_fasta(open(self.multifasta))
		snps_ordered_nonhom, id_snps_nonhom = findSNPs.setupSNPdict(nonhomoplasious_sites, alignment)
		with open(str(self.output_prefix)+'.nonhomoplasious.fasta', 'w') as fp:
			for i in id_snps_nonhom:
				if id_snps_nonhom[i].count('-') > round(len(id_snps_nonhom[i])*0.45):
					self.logger.warning(str(i)+' not in non homoplasious fasta!')
				else:
					fp.write('>'+ i +'\n')
					fp.write(''.join(id_snps_nonhom[i])+'\n')
		subprocess.call([self.fasttree_exec+' -nt '+str(self.output_prefix)+'.nonhomoplasious.fasta > '+str(self.output_prefix)+'.nonhomoplasious.nwk'], stdout=subprocess.DEVNULL, shell=True)
		
		peaks = set(positions) - set(filtered_positions)
		with open(str(self.output_prefix) + '_homoplasious_sites.txt', 'w') as f:
			for i in homoplasious_sites:
				f.write(str(i) +'\n')
		GenesWithclusteredSNPs = self.findFeatures(gff, homoplasious_sites)
		try:
			CodingGenes = SeqIO.parse(self.faa_file, 'fasta')
			featuredict = self.featuresdict(gff)
			PhageGenes = []
			for gene in CodingGenes:
				if str(gene.id) in GenesWithclusteredSNPs:
					self.logger.warning(gene.id)
					with open("temp.fasta", "w") as output_handle:
						SeqIO.write(gene, output_handle, 'fasta')
						
					blastp_cline = NcbiblastpCommandline(query="temp.fasta", db=pkg_resources.resource_filename(__name__, 'databases/prophage_virus.db'), evalue=1.1E-11,
														 outfmt=5, out="test.xml", max_hsps=1, max_target_seqs=1, num_threads=4)
					blastp_cline()
					alignments = SearchIO.parse('test.xml', 'blast-xml')
					for alignment in alignments:
						for hit in alignment.hits:
							if hit.id.startswith('PHAGE') or hit.id.startswith('PROPHAGE'):
								PhageGenes.append(gene.id)
								self.logger.warning(hit.id[0:75] + '...')
			phage_ranges = [featuredict[x] for x in PhageGenes]
			phagepositions = []
			for p in phage_ranges:
				phagepositions += list(range(p[0], p[1] + 1))
			variable_phage_positions = sorted(list(set(positions) & set(phagepositions)))
			with open(str(self.output_prefix) + '_phage_SNPs.txt', 'w') as f:
				for i in variable_phage_positions:
					f.write(str(i) + '\n')
			window_phage, window_start = self.sliding_window(n, variable_phage_positions)
			sns.plt.plot(window_start, window_phage, 'go')
			
		except FileNotFoundError:
			self.logger.warning('Can not find AA gene file; will not perform phage related gene scan. Please provide .faa file if you\n want to perform this procedure. ')
		
		sns.plt.plot(window_start, window_hom, 'ro')
		sns.plt.plot(window_start, window_count, color='b', alpha=0.7)
		plt.savefig(str(self.output_prefix) + 'SNP_homoplasy_density.pdf')
