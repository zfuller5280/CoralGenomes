from pyfaidx import Fasta
import sys
import random
import subprocess
from subprocess import Popen, PIPE

ref_seq = sys.argv[1]
input_vcf = sys.argv[2]

ref_seq = Fasta(ref_seq)

variants = []
with open(input_vcf) as fp:
    for line in fp:
        if line.startswith("#"):
            if line.startswith("##bcftools_viewCommand=view -r"):
                line = line.split(" ")
                region = line[2]
                chrom, start, stop = region.split(":")[0], int(region.split(":")[1].split("-")[0]), int(region.split(":")[1].split("-")[1])
            elif line.startswith("#CHROM"):
                header =  line.strip("\r\n").split("\t")
        else:
            line = line.strip("\r\n").split("\t")
            variants.append(line)

region_seq = ref_seq[chrom][start-1:stop]
samps = header[9:]

hap_choices = []
for _ in range(len(samps)):
    k = random.randint(0, 1)
    hap_choices.append(k)

hap_seqs = [[] for i in range(len(samps))]

vcf_count = 0
num_vars = len(variants)
for n, site in enumerate(region_seq):
    site = str(site)
    if vcf_count >= num_vars:
        for x in hap_seqs:
            x.append(str(site))
    else:
        if n == (int(variants[vcf_count][1]) - start):
            ref, alt = variants[vcf_count][3], variants[vcf_count][4]
            #print site, ref, alt
            for ind, mut in enumerate(variants[vcf_count][9:]):

                geno = mut.split("|")
                if geno[hap_choices[ind]] == "0":
                    hap_seqs[ind].append(str(ref))
                else:
                     hap_seqs[ind].append(str(alt))
            vcf_count += 1
        else:

            for x in hap_seqs:

                x.append(str(site))

def get_outgroup_seq(outgroup, chrom, start, stop):
    cmd1 = ["/Users/zachfuller/anaconda/bin/samtools","faidx","/Users/zachfuller/Amil_v2.01/Amil.v2.01.chrs.fasta","%s:%i-%i"%(chrom, start, stop)]
    cmd2 = ["/usr/local/ncbi/blast/bin/blastn","-query","-","-db","/Users/zachfuller/%s"%(outgroup),"-outfmt","6 sseq"]
    p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=subprocess.PIPE)
    output = p2.stdout.readlines()
    return output[0].strip("\n")

outgroups = ["Adigit","Aten" d]
for n, outgroup in enumerate(outgroups):
    outgroup_seq = get_outgroup_seq(outgroup, chrom, start, stop)
    print ">%s"%outgroups[n]
    print outgroup_seq

for n, x in enumerate(hap_seqs):
    print ">%s"%samps[n]
    print ''.join(x)
