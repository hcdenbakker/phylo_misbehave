# phylo_misbehave
Pipeline to infer SNP-dense regions and the influence of different filtering strategies. Work in progress, hence the bare bones description! I am just sharing this in case it may be useful to anyone.




# Installation
Requirements:

-Python 3 

-Ete3 (http://etetoolkit.org/); not need to install ete3_external_apps

-Biopython

-FastTree (http://www.microbesonline.org/fasttree/), remember to compile get the version -DUSE_DOUBLE (http://darlinglab.org/blog/2015/03/23/not-so-fast-fasttree.html)

-NCBI blast


## Debian Testing/ Ubuntu 16.04
These instructions assume you have root permissions:

```
sudo apt-get update -qq
sudo apt-get install -y git python3 python3-setuptools python3-biopython python3-pip ncbi-blast+ gcc cython3 python3-dev fasttree
pip3 install git+git://github.com/hcdenbakker/phylo_misbehave.git
```

## Docker
Install [Docker](https://www.docker.com/).  We have a docker container which gets automatically built from the latest version of PhyloMisbehave.  To use it you would use a command such as this (substituting in your directories), where your files are assumed to be stored in /home/ubuntu/data:
```
docker run --rm -it -v /home/ubuntu/data:/data hcdenbakker/phylo_misbehave your_multi_fasta.fasta prefix_for_output_files gff_file.gff faa_file.faa
```

# How to use:
```
usage: phylomisbehave [options] your_multi_fasta.fasta prefix_for_output_files gff_file.gff faa_file.faa

Infer SNP-dense regions and the influence of different filtering strategies

positional arguments:
  multifasta            MultiFASTA file
  output_prefix         Output prefix
  gff_file              GFF file
  faa_file              FAA file corresponding to GFF file

optional arguments:
  -h, --help            show this help message and exit
  --threads THREADS, -t THREADS
                        Number of threads (default: 1)
  --evalue EVALUE, -e EVALUE
                        Blast E-value (default: 1.1e-11)
  --min_snps MIN_SNPS, -m MIN_SNPS
                        Min SNPs in a sliding window to consider (default: 5)
  --verbose, -v         Turn on more debugging output (default: False)
  --version             show program's version number and exit
```
gff/faa files are used to do a scan of prophage/phage related genes in the reference genome. Both can be generated by Prokka (https://github.com/tseemann/prokka).

### output:

variable sites only multifasta (created with a partly rewitten version of https://github.com/bewt85/PySnpSites)

tree based on unfiltered multifasta

variable sites multifasta with snp dense regions filtered

tree based on snp dense regions filtered dataset

multifasta without homoplasious sites

tree based on homoplasious site filtered dataset

multifasta with prophage/phage related snps filtered + tree (using prophage and virus database found here: http://phast.wishartlab.com/Download.html)

list of SNPs in phage-related regions

list of SNPs inferred to be homoplasious
