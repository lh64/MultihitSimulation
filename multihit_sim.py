import numpy as np
import matplotlib.pyplot as plt
import re

filename = 'INPUT/gene_length_mut_5strains.csv'
suffix = re.match(".+mut(.*)\.csv",filename).group(1)
data = np.genfromtxt(filename, dtype=None, delimiter=',', names=True, encoding='UTF-8')

print(data)
print(data.dtype.names)

HO_nonsyn = [(d['Length'], 2*d['Nonsyn']) for d in data 
             if d['Orientation'] == 'HO' and d['Length'] > 200]
CD_nonsyn = [(d['Length'], d['Nonsyn']) for d in data 
             if d['Orientation'] == 'CD' and d['Length'] > 200]

np.random.seed(0)

niter = 10000

sum_hsmut_sites_HO = [0]*niter
sum_hsmut_sites_CD = [0]*niter
sum_hsmut_genes_HO = [0]*niter
sum_hsmut_genes_CD = [0]*niter
sum_parmut_HO  = [0]*niter
sum_parmut_CD  = [0]*niter

for n in range(niter):
    print(n)
    sum_hsmut_sites = [0,0] # [HO,CD]
    sum_hsmut_genes = [0,0]
    sum_parmut = [0,0]
    for i,type in enumerate([HO_nonsyn, CD_nonsyn]):
        for gene in type:
            x = np.random.randint(gene[0], size=gene[1])
            pos, hits = np.unique(x, return_counts=True)
            '''
            while (len(pos) < gene[1]):
                x = np.append(x, np.random.randint(gene[0], size=gene[1]-len(pos)))
                pos, hits = np.unique(x, return_counts=True)
            '''
            multihits = [h for h in hits if h > 1]
            n_multihits = len(multihits)

            if n_multihits > 0:
                sum_hsmut_sites[i] += n_multihits
                sum_hsmut_genes[i] += 1
                # Check for parallel mutations
                for mh in multihits:
                    # 20 amino acids to choose from (there are 21 total)
                    xx = np.random.randint(20, size=mh) 
                    aa, nmut = np.unique(xx, return_counts=True)
                    parallel = [nm for nm in nmut if nm > 1]
                    n_parallel = len(parallel)
                    if n_parallel > 0:
                        sum_parmut[i] += n_parallel
            
    sum_hsmut_sites_HO[n] = sum_hsmut_sites[0]
    sum_hsmut_sites_CD[n] = sum_hsmut_sites[1]
    sum_hsmut_genes_HO[n] = sum_hsmut_genes[0]
    sum_hsmut_genes_CD[n] = sum_hsmut_genes[1]
    sum_parmut_HO[n]  = sum_parmut[0]
    sum_parmut_CD[n]  = sum_parmut[1]

with open('OUTPUT/hsmut_sites%s.csv' % suffix, 'w') as f:
    f.write("HO,CD\n")
    for n in range(niter):
        f.write("%d,%d\n" % (sum_hsmut_sites_HO[n],sum_hsmut_sites_CD[n]))

with open('OUTPUT/hsmut_genes%s.csv' % suffix, 'w') as f:
    f.write("HO,CD\n")
    for n in range(niter):
        f.write("%d,%d\n" % (sum_hsmut_genes_HO[n],sum_hsmut_genes_CD[n]))

with open('OUTPUT/parmut%s.csv' % suffix, 'w') as f:
    f.write("HO,CD\n")
    for n in range(niter):
        f.write("%d,%d\n" % (sum_parmut_HO[n],sum_parmut_CD[n]))
