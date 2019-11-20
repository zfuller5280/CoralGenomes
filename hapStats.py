from optparse import OptionParser
from collections import Counter
import itertools

USAGE= """Usage: %prog [options] -i infile.txt"""
OPT_DEFAULTS={'infile':'-'}
DESCRIPTION="""Program description"""
EPILOG="""Requirements:"""

def get_options(defaults, usage, description='',epilog=''):
    """Get options, print usage text."""
    parser = OptionParser(usage=usage,description=description,epilog=epilog)
    parser.add_option("-i","--infile",action="store",dest="infile",type="string",
                      default=defaults.get('infile'),
                      help='Name of input file')
    parser.add_option("-w","--window_size",action="store",dest="window_size",type="int",
                      default=defaults.get('window'),
                      help='Size of window')

    (options,args)=parser.parse_args()

    return (options, args)

def read_vcf_windows(vcf, window_size):
    if window_size % 2 == 0:
        window_size += 1

    middleIndex = (window_size - 1)/2
    window, pos_L = [], []
    start, stop = 1, window_size
    current_chrom = None
    with open(vcf) as fp:
        for line in fp:
            if line.startswith("#"):
                continue
            line = line.strip("\r\n").split("\t")
            chrom, pos = line[0], int(line[1])

            if current_chrom != chrom:
                if len(window) > 0: yield current_chrom, window, start, stop, focal_pos
                window, pos_L = [], []
                start = pos
                current_chrom = chrom
                window.append(line)
                pos_L.append(pos)
                continue
            if len(window) > window_size:
                stop = pos
                focal_pos = pos_L[middleIndex]
                yield current_chrom, window, start, stop, focal_pos
                window.pop(0)
                pos_L.pop(0)
                window.append(line)
                pos_L.append(pos)
                start = pos_L[0]
            else:
                window.append(line)
                pos_L.append(pos)

def calc_hap_stats_in_window(window):
    window_snps = []
    for site in window:
        ref, alt = site[3:5]
        snps = site[9:]
        snps = [w.replace('.', ref) for w in snps]
        snps = [int(y) for x in snps for y in x.split("|")]
        snps = [ref if x==0 else alt for x in snps]
        window_snps.append(snps)
    sites = map(list, zip(*window_snps))
    haps = [''.join(x) for x in sites]
    hap_counts = Counter(haps).most_common()
    count_h1, count_h2 = hap_counts[0][1],hap_counts[1][1]
    total_uniq_haps = len(hap_counts)
    total_inds = len(sites)
    p1, p2 = float(count_h1)/total_inds, float(count_h2)/total_inds
    h1 = sum([(float(count_h[1])/total_inds)**2 for count_h in hap_counts])
    h12 = h1 + (2*p1*p2)
    h2 = h1 - (p1**2)
    h2_1 = h2/h1
    freq_spec = [x[1] for x in hap_counts]
    hap_pi = float(sum([x[0] * x[1] for x in list(itertools.combinations(freq_spec, 2))]))/float((total_inds*(total_inds-1))/2)
    weight_hap_pi = hap_pi * (float(1)/(total_uniq_haps-1))
    return total_inds, total_uniq_haps, p1, p2, h1, h12, h2, h2_1, hap_pi, weight_hap_pi


def main():
    (options,args) = get_options(OPT_DEFAULTS, USAGE, DESCRIPTION, EPILOG)
    infile = options.infile
    window_size = options.window_size

    for chrom, window, start, stop, focal_pos in read_vcf_windows(infile, window_size):
        total_inds, total_uniq_haps, p1, p2, h1, h12, h2, h2_1, hap_pi, weight_hap_pi  = calc_hap_stats_in_window(window)
        print chrom, start, stop, focal_pos, total_inds, total_uniq_haps, p1, p2, h1, h12, h2, h2_1, hap_pi, weight_hap_pi
        #break
        #break

if __name__ == '__main__':
    main()
