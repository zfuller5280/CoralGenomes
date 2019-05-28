import numpy as np
import pandas as pd
import sys


vcf_path = sys.argv[2]
truth_set = sys.argv[1]

sample_name = sys.argv[3]

sample_vcf = pd.read_table(vcf_path, comment="#", header=None)

def get_sample_idx(truth_set, sample_name):
    with open(truth_set) as fp:
        for line in fp:
            if line.startswith("#CHROM"):
                line = line.strip("\r\n").split("\t")
                sample_idx = line.index(sample_name)
                break
    return sample_idx

sample_vcf.rename(columns={sample_vcf.columns[4]:'alt'}, inplace=True)
sample_vcf.rename(columns={sample_vcf.columns[1]:'pos_raw'}, inplace=True)
sample_vcf.rename(columns={sample_vcf.columns[0]:'chrom_raw'}, inplace=True)
sample_vcf.rename(columns={sample_vcf.columns[7]:'info'}, inplace=True)
sample_vcf["gt"] = (sample_vcf.ix[:,9].str.split(":").str.get(0))
sample_vcf['chrom']=sample_vcf['chrom_raw'].str.split(':', 1).str[0]
sample_vcf['pos_str']=sample_vcf['chrom_raw'].str.split(':', 1).str[1]
sample_vcf['pos']=sample_vcf['pos_str'].str.split('-',1).str[0].astype(np.int64)
sample_vcf['depth']=sample_vcf['info'].str.split(';',1).str[0].str.split('=',1).str[1].astype(np.int64)
sample_vcf['pos']=sample_vcf['pos'] + (sample_vcf['pos_raw']-1)

sample_vcf = sample_vcf[sample_vcf['depth']>=50]


sample_idx = get_sample_idx(truth_set, sample_name)
#print sample_vcf
alt_mismatch = 0
current_chrom = None
gt_dict = {}
window_length = 0
with open(truth_set) as fp:
    for line in fp:
        if line.startswith("#"): continue
        line = line.strip("\r\n").split("\t")
        sample_gt = line[sample_idx].split(":")[0]
        chrom, pos = line[0:2]
        pos = int(pos)
        if chrom != current_chrom:
            current_chrom = chrom
            window_length = 0
            chrom_vcf = sample_vcf.ix[(sample_vcf["chrom"])==current_chrom]
            #chrom_vcf = sample_vcf.ix[(sample_vcf["chrom"])==current_chrom]

        site = chrom_vcf.ix[(chrom_vcf["chrom"]==chrom) & (chrom_vcf["pos"]==int(pos))]
        if site.empty:
            site_gt = "./."
        else:
            site_gt = site["gt"].values[0]
            if "1" in site_gt:

                alt = site["alt"].str.split(",").str.get(0).values[0]

                # if alt != line[4]:
                #     alt_mismatch += 1
                #     #print alt, line[4]
                #     continue
        #print site
        gt_key = "%s_%s"%(sample_gt, site_gt)
        if not gt_key in gt_dict:
            gt_dict[gt_key] = 1
        else:
            gt_dict[gt_key] += 1
        if site_gt != "./." and sample_gt != "./.":
            #print line
            print sample_name, sample_gt, site_gt, pos, line[5], window_length, line[7], line[sample_idx]
gt_dict["MISMATCH"]=alt_mismatch
gt_dict["SAMPLE"]=sample_name
keys = gt_dict.keys()
vals = gt_dict.values()
#print "\t".join(keys)
#print "\t".join([str(c) for c in vals])
        #break
