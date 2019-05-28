import numpy as np
import pandas as pd
import sys


vcf_path = sys.argv[2]
truth_set = sys.argv[1]

#sample_name = sys.argv[3]



def get_vcf_header(truth_set):
    with open(truth_set) as fp:
        for line in fp:
            if line.startswith("#CHROM"):
                header = line.strip("\r\n").split("\t")
                break
    return header

header = get_vcf_header(truth_set)
header = [s.replace('.sorted.bam', '').split("_")[0] for s in header]

truth_vcf = pd.read_table(truth_set, comment="#", header=None)

def get_sample_idx(truth_set, sample_name):
    with open(truth_set) as fp:
        for line in fp:
            if line.startswith("#CHROM"):
                line = line.strip("\r\n").split("\t")
                sample_idx = line.index(sample_name)
                break
    return sample_idx

truth_vcf.columns = header
truth_vcf['chrom']=truth_vcf['#CHROM'].str.split(':', 1).str[0]
truth_vcf['pos_str']=truth_vcf['#CHROM'].str.split(':', 1).str[1]

truth_vcf['pos']=truth_vcf['pos_str'].str.split('-',1).str[0].astype(np.int64)
truth_vcf['pos']=truth_vcf['pos'] + (truth_vcf['POS']-1)

sample_header = get_vcf_header(vcf_path)

site_stat = []
with open(vcf_path) as fp:
    for line in fp:
        if line.startswith("#"): continue
        line = line.strip("\r\n").split("\t")
        chrom, pos = line[0:2]
        pos = int(pos)
        #print pos
        site = truth_vcf.ix[truth_vcf["pos"]==int(pos)]
        if site.empty:
            continue
        sample_gts, truth_gts = [], []
        for sample in sample_header[9:]:
            sample_gts.append(line[sample_header.index(sample)].split(":")[0])
            truth_gts.append(str(site[sample].values[0]).split(":")[0])
        nonmissing_truth = [x for x in truth_gts if x != "./."]
        #print nonmissing_truth
        nonmissing_truth = [int(x) for xs in nonmissing_truth for x in xs.split('/')]
        nonmissing_truth_alts = sum(nonmissing_truth)
        if nonmissing_truth_alts < 1:
            continue
        sample_missing_idx = [i for i, x in enumerate(sample_gts) if x == "./."]
        sample_gts = [i for j, i in enumerate(sample_gts) if j not in sample_missing_idx]
        truth_gts = [i for j, i in enumerate(truth_gts) if j not in sample_missing_idx]

        truth_missing_idx = [i for i, x in enumerate(truth_gts) if x == "./."]
        truth_gts = [i for j, i in enumerate(truth_gts) if j not in truth_missing_idx]
        sample_gts = [i for j, i in enumerate(sample_gts) if j not in truth_missing_idx]

        samp_freqs = [int(x) for xs in sample_gts for x in xs.split('/')]
        truth_freqs = [int(x) for xs in truth_gts for x in xs.split('/')]
        samp_freq, truth_freq = float(sum(samp_freqs))/len(samp_freqs), float(sum(truth_freqs))/len(truth_freqs)

        c = (np.array(truth_gts) == np.array(sample_gts))
        if False in c:
            stat = 0
        else:
            stat = 1
        site_stat.append(stat)
        print stat, samp_freq, truth_freq, line[7]
#print float(sum(site_stat))/len(site_stat)
