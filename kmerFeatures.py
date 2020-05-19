#from optparse import OptionParser
from collections import Counter
from kpal.klib import Profile
import numpy as np
import pandas as pd
import Queue as queue
import multiprocessing as mp

USAGE= """Usage: %prog [options] -i infile.txt"""
OPT_DEFAULTS={'infile':'-','test_file':'-'}
DESCRIPTION="""Program description"""
EPILOG="""Requirements:"""

def get_options(defaults, usage, description='',epilog=''):
    """Get options, print usage text."""
    parser=OptionParser(usage=usage,description=description,epilog=epilog)
    parser.add_option("-i","--infile",action="store",dest="infile",type="string",
                      default=defaults.get('infiles'),
                      help='Name of input file of contigs, in .fasta')
    parser.add_option("-k","--kmers",action="store",dest="kmers",type="string",
                      default=defaults.get('kmers'),
                      help='Sizes of k-mers to use as features, comma separated list')
    (options,args)=parser.parse_args()

    return (options, args)

def unique(seq, idfun=None):
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]

##FASTA reader function
def read_fasta(fp):
    name, seq=None,[]
    for line in fp:
        line=line.rstrip()
        if line.startswith(">"):
            if name: yield(name,''.join(seq))
            name,seq=line,[]
        else:
            seq.append(line)
    if name: yield(name,''.join(seq))

def kmer_count(seq, k):
    f = {}
    for x in range(len(seq) + 1-k):
        kmer = seq[x:x+k]
        f[kmer] = f.get(kmer, 0) + 1
    return f

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])

def build_feature_matrix(kmer_feature_names, kmer_features_matrices, names):
    kmer_dfs = []
    for n, i in enumerate(kmer_feature_names):
        df = pd.DataFrame(data = kmer_features_matrices[n], columns = kmer_feature_names[n], index = names)
        kmer_dfs.append(df)

    kmer_df = reduce(lambda x, y: pd.merge(x, y, left_index=True, right_index=True, how='outer'), kmer_dfs)
    return kmer_df

def kmer_feature_matrix(names, seqs, kmer_sizes, q):
    freq_L, keys_L = [[] for i in range(len(kmer_sizes))], [[] for i in range(len(kmer_sizes))]

    for name, seq in zip(names, seqs):
        seq += reverse_complement(seq)
        for n, k in enumerate(kmer_sizes):
            p = Profile.from_sequences([seq], int(k))

            #print p.balance
            f = {}
            for i in range(p.number):
                kmer = p.binary_to_dna(i)
                f[kmer] = p.counts[i]
            kmer_dict = f

            freqs = [float(x)/sum(kmer_dict.values()) for x in kmer_dict.values()]
            keys = kmer_dict.keys()
            freq_L[n].append(freqs), keys_L[n].append(keys)
        print name
    kmer_feature_names = []
    kmer_features_matrices = []

    for n, k in enumerate(kmer_sizes):
        flatten = lambda l: [item for sublist in l for item in sublist]
        kmers = unique(flatten(keys_L[n]))

        kmer_features = np.zeros(shape=(len(names),len(kmers)))
        for seq_n, seq in enumerate(freq_L[n]):
            for kmer_n, kmer in enumerate(keys_L[n][seq_n]):
                kmer_features[seq_n][kmers.index(kmer)] = seq[kmer_n]
        kmer_features_matrices.append(kmer_features)
        kmer_feature_names.append(kmers)
    kmer_feature_matrix = build_feature_matrix(kmer_feature_names, kmer_features_matrices, names)
    q.put(kmer_feature_matrix)
