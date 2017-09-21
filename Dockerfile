FROM debian:testing
MAINTAINER ap13@sanger.ac.uk

RUN apt-get update -qq && apt-get install -y git python3 python3-setuptools python3-biopython python3-pip ncbi-blast+ gcc cython python3-dev fasttree

RUN pip3 install git+git://github.com/andrewjpage/phylo_misbehave.git

WORKDIR /data