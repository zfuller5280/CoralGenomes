import sys
import subprocess
from subprocess import Popen, PIPE

def read_intervals(intervals):
    interval_L = []
    with open(intervals) as fp:
        for line in fp:
            line = line.strip("\r\n").split("\t")
            chrom, start, stop = line[0], int(line[1]), int(line[2])
            interval_L.append([chrom, start, stop])
    return interval_L

def parse_site_pi(site_pi):
    pi = 0
    for i in site_pi[1:]:
        site = i.strip("\r\n").split("\t")
        pi += float(site[-1])
        #print site
    return pi

def get_uncallable_sites(chrom, start, stop):
    cmd1 = ["bcftools","view","-r","%s:%i-%i"%(chrom, start, stop),"Millepora_all_GT_filt.vcf.gz"]
    cmd2 = ["bcftools","query","-f",""'%POS\n'",-"]
    p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE)
    callable_sites = len(p2.stdout.readlines())
    uncallable_sites = (stop-start) - callable_sites
    return uncallable_sites, callable_sites


def main():
    input_vcf = sys.argv[1]
    intervals = sys.argv[2]

    interval_L = read_intervals(intervals)
    for interval in interval_L:
        chrom, start, stop = interval
        cmd1 = ["bcftools","view","-r","%s:%i-%i"%(chrom, start, stop),input_vcf]
        cmd2 = ["vcftools","--vcf","-","--site-pi","--stdout"]
        p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE)
        site_pi = p2.stdout.readlines()
        uncallable_sites, callable_sites = get_uncallable_sites(chrom, start, stop)
        print chrom, start, stop, parse_site_pi(site_pi)/((stop-start)-uncallable_sites), callable_sites
        #break



if __name__ == '__main__':
    main()
