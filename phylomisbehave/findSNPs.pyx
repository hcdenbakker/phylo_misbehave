from collections import OrderedDict

def parse_fasta(input_fasta):
    '''fast fasta parser from Ben Taylor
        https://github.com/bewt85/PySnpSites/blob/master/snp_sites_extensions.pyx'''
    for line in input_fasta:
        if line[0] == '>':
            break
    sequence_name = line[1:].rstrip()
    sequence_lines = []
    for line in input_fasta:
        if line[0] == '>':
            yield (sequence_name, "".join(sequence_lines))
            sequence_name = line[1:].rstrip()
            sequence_lines = []
        else:
            sequence_lines.append(line.rstrip())
    yield(sequence_name, "".join(sequence_lines))

def bestref(alignment):
    min_amb = 15000000
    cdef str record, seq
    for record, seq in alignment:
        if seq.count('-') < min_amb:
            min_amb = seq.count('-')
            best_ref = record
            best_seq = seq
    return (best_ref, best_seq)

def findSNPs(refseq, alignment):
    varsites = []
    cdef int i
    cdef str r, s
    for record, seq in alignment:
        #print(record)
        for i in range(len(refseq)):
            r, s  = refseq[i], seq[i]
            if r != s:
                varsites.append(i)
    return sorted(list(set(varsites)))

def setupSNPdict(VariableSites, alignment):
    snps_raw = {}
    id_snps ={}
    cdef int i, index
    cdef str id, s, seq
    cdef tuple tup
    sequences =[]
    ids = []
    for id, s in alignment:
        sequences.append(s)
        ids.append(id)
        id_snps[id] = []
    for i in VariableSites:
        characters =[]
        id_char =[]
        for index, seq in enumerate(sequences):
            characters.append(seq[i])
            id_char.append((ids[index], seq[i]))
        #print(characters)
        if len(set(characters) - set('-') - set('N') - set('n')) > 1:
            snps_raw[i] = id_char
            for tup in id_char:
                id_snps[tup[0]].append(tup[1])
    snps_ordered = OrderedDict(sorted(snps_raw.items(), key=lambda t: t[0]))
    return snps_ordered, id_snps