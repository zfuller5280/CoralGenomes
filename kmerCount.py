from optparse import OptionParser
from kmerFeatures import *
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
    parser.add_option("-p","--processors",action="store",dest="processors",type="int",
                      default=defaults.get('processors'),
                      help='Number of processors to use')
    parser.add_option("-o","--outfile",action="store",dest="outfile",type="string",
                      default=defaults.get('outfile'),
                      help='Name of file to write output .csv to')

    (options,args)=parser.parse_args()

    return (options, args)

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

def chunk(xs, n):
    '''Split the list, xs, into n chunks'''
    L = len(xs)
    assert 0 < n <= L
    s, r = divmod(L, n)
    chunks = [xs[p:p+s] for p in range(0, L, s)]
    chunks[n-1:] = [xs[-r-s:]]
    return chunks

def main():
    (options,args)=get_options(OPT_DEFAULTS, USAGE, DESCRIPTION, EPILOG)
    infile = options.infile
    outfile = options.outfile
    kmer_sizes = options.kmers.split(",")
    processors = options.processors

    print "Reading input sequences..."
    names, seqs = [], []
    with open(infile) as fp:
        for name, seq in read_fasta(fp):
            seq = seq.upper()

            names.append(name), seqs.append(seq)

    ##Setup multiprocessing for kmer frequency distribution
    chunked_names, chunked_seqs = chunk(names, processors), chunk(seqs, processors)
    t_c = 0
    processes_L=[]

    queues=[mp.Queue() for i in chunked_names]

    for n, i in enumerate(chunked_names):
        processes_L.append(mp.Process(target=kmer_feature_matrix, args=(i, chunked_seqs[n], kmer_sizes, queues[t_c])))
        t_c += 1

    process_results=[]

    #Start the processes in the job queue
    for p in processes_L:
        p.start()

    #Return the output from each process in the queue and place in list
    for q in queues:
        process_results.append(q.get())

    #Coordinate and join each process
    for p in processes_L:
        p.join()

    kmer_features = pd.concat(process_results)
    kmer_features['GC_content'] = kmer_features['G'] + kmer_features['C']
    print kmer_features
    kmer_features.to_csv(outfile)

if __name__ == '__main__':
    main()
