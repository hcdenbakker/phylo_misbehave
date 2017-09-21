# phylo_misbehave
Pipeline to infer SNP-dense regions and the influence of different filtering strategies. Work in progress, hence the bare bones description! I am just sharing this in case it may be useful to anyone.

Requirements:

-Python 3 

-Ete3 (http://etetoolkit.org/); not need to install ete3_external_apps

-Biopython

-FastTree (http://www.microbesonline.org/fasttree/), remember to compile get the version -DUSE_DOUBLE (http://darlinglab.org/blog/2015/03/23/not-so-fast-fasttree.html)

-NCBI blast

Probably forgot some requirements, will add later!

##How to use:

python misbehave.py your_multi_fasta.fasta prefix_for_output_files gff_file.gff

gff file is used to do a scan of prophage/phage related genes in the reference genome and assumes a .faa file of the genes in genome is present in the execution folder. Both can be generated by Prokka (https://github.com/tseemann/prokka).

##output:

variable sites only multifasta (created with a partly rewitten version of https://github.com/bewt85/PySnpSites)

tree based on unfiltered multifasta

variable sites multifasta with snp dense regions filtered

tree based on snp dense regions filtered dataset

multifasta without homoplasious sites

tree based on homoplasious site filtered dataset

multifasta with prophage/phage related snps filtered + tree (using prophage and virus database found here: http://phast.wishartlab.com/Download.html)

list of SNPs in phage-related regions
list of SNPs inferred to be homoplasious
