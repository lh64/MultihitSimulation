import numpy as np
import matplotlib.pyplot as plt
import re
from cmath import nan

########## READ SIMULATION OUTPUT FILE ##########
def make_plots(SUFFIX,EC,FC,LW,LABEL):
    
    for suffix,ec,fc,lw,label in zip(SUFFIX,EC,FC,LW,LABEL):
        
        data_exp = np.genfromtxt('INPUT/gene_length_mut%s.csv' % suffix, 
                             dtype=None, delimiter=',', names=True, encoding='UTF-8')
        
        HO_mut_exp = [(d['HSmut'], d['Convergent'], d['Parallel'], d['Coincidental']) 
                  for d in data_exp if d['Orientation'] == 'HO' and d['Length'] > 200]
        CD_mut_exp = [(d['HSmut'], d['Convergent'], d['Parallel'], d['Coincidental']) 
                  for d in data_exp if d['Orientation'] == 'CD' and d['Length'] > 200]
        
        n_ho_genes = len(HO_mut_exp)
        n_cd_genes = len(CD_mut_exp)
                
        ho_hsmut_exp = [x[0] for x in HO_mut_exp if x[0]>0]
        cd_hsmut_exp = [x[0] for x in CD_mut_exp if x[0]>0]
        
        print(suffix)
        print(sum(ho_hsmut_exp), n_ho_genes, '%.2f' % (sum(ho_hsmut_exp)/n_ho_genes))
        print(sum(cd_hsmut_exp), n_cd_genes, '%.2f' % (sum(cd_hsmut_exp)/n_cd_genes))
        
        ratio_sites_exp = sum(ho_hsmut_exp)/sum(cd_hsmut_exp)
        ratio_genes_exp = len(ho_hsmut_exp)/len(cd_hsmut_exp)
        
#         suffix += '_test'
        
        data_sim_sites = np.genfromtxt('sum_hsmut_sites%s.csv' % suffix, 
                            dtype=None, delimiter=',', names=True, encoding='UTF-8')
        data_sim_genes = np.genfromtxt('sum_hsmut_genes%s.csv' % suffix, 
                            dtype=None, delimiter=',', names=True, encoding='UTF-8')
        
        ho_hsmut_sites_sim = [d[0] for d in data_sim_sites]
        cd_hsmut_sites_sim = [d[1] for d in data_sim_sites]
        ratio_sites_sim = [d[0]/d[1] for d in data_sim_sites]
        
        ho_hsmut_genes_sim = [d[0] for d in data_sim_genes]
        cd_hsmut_genes_sim = [d[1] for d in data_sim_genes]
        ratio_genes_sim = [d[0]/d[1] for d in data_sim_genes]


        
        plt.figure('R_N')
        dx = 0.025 #0.1
#         bins = np.arange(0, 2.5+dx, dx)
        bins = np.arange(0, 0.8+dx, dx)
        plt.hist(ratio_sites_sim, bins=bins, density=False, 
                 weights=[1/len(ratio_sites_sim)]*len(ratio_sites_sim), 
                 color='white', ec=ec, fc=fc, lw=lw, label=label)
        plt.plot([ratio_sites_exp], [0.2], '*', ms=12, mec=fc, mfc=ec, label='observed')
        plt.xlabel(r'$R_N$')
        plt.ylabel('frequency')
#         plt.xticks([bins[i] for i in range(len(bins)) if np.mod(i,4)==0])
        plt.ylim(ymax=0.5) #0.5
        plt.legend(loc=0)

        plt.figure('R_G')
        dx = 0.025 #0.1
#         bins = np.arange(0, 2+dx, dx)
        bins = np.arange(0, 0.8+dx, dx)
        plt.hist(ratio_genes_sim, bins=bins, density=False, 
                 weights=[1/len(ratio_genes_sim)]*len(ratio_genes_sim), 
                 color='white', ec=ec, fc=fc, lw=lw, label=label)
        plt.plot([ratio_genes_exp], [0.2], '*', ms=12, mec=fc, mfc=ec, label='observed')
        plt.xlabel(r'$R_G$')
        plt.ylabel('frequency')
#         plt.xticks([bins[i] for i in range(len(bins)) if np.mod(i,4)==0])
        plt.ylim(ymax=0.5) #0.5
        plt.legend(loc=0)
        
        def get_bins(x, y):
            dmin = min(min(x),min(y))
            dmax = max(max(x),max(y))
            order = 0
            n = 2
            while dmax > 2*10**(order+(order+1)*n):
                order += 1
            step = 10**order
            if dmin % 10*step != 0:
                dmin = int(round((dmin-5*step)/(10*step))*(10*step))
            if dmax % 10*step != 0:
                dmax = int(round((dmax+5*step)/(10*step))*(10*step))
            bins = range(dmin,dmax,step)
            return bins
        
        bins = get_bins(ho_hsmut_sites_sim, cd_hsmut_sites_sim)
        
        plt.figure('hsmut_sites%s' % suffix)
        plt.hist(ho_hsmut_sites_sim, bins=bins, color='white', fc='0.75', 
                 ec='red', lw=1, label='simulated HO')
        plt.hist(cd_hsmut_sites_sim, bins=bins, color='white', fc='None', 
                 ec='blue', lw=2, label='simulated CD')
        plt.annotate('Observed HO = {:,}'.format(sum(ho_hsmut_exp)) , xy=(0.6,0.75), 
                     xycoords='axes fraction', weight='bold', color='red')
        plt.annotate('Observed CD = {:,}'.format(sum(cd_hsmut_exp)), xy=(0.6,0.69), 
                     xycoords='axes fraction', weight='bold', color='blue')
        plt.xlabel('# multi-hit sites') 
        plt.ylabel('number')
        # plt.xticks(range(0,11))
        plt.legend(loc=0)        
        
        bins = get_bins(ho_hsmut_genes_sim, cd_hsmut_genes_sim)
        
        plt.figure('hsmut_genes%s' % suffix)
        plt.hist(ho_hsmut_genes_sim, bins=bins, color='white', fc='0.75', 
                 ec='red', lw=1, label='simulated HO')
        plt.hist(cd_hsmut_genes_sim, bins=bins, color='white', fc='None', 
                 ec='blue', lw=2, label='simulated CD')
        plt.annotate('Observed HO = {:,}'.format(len(ho_hsmut_exp)) , xy=(0.6,0.75), 
                     xycoords='axes fraction', weight='bold', color='red')
        plt.annotate('Observed CD = {:,}'.format(len(cd_hsmut_exp)), xy=(0.6,0.69), 
                     xycoords='axes fraction', weight='bold', color='blue')
        plt.xlabel('# multi-hit genes') 
        plt.ylabel('number')
        # plt.xticks(range(0,11))
        plt.legend(loc=0)    
        
        ##### PARALLEL MUTATIONS #####
        
        ho_par_exp = [x[2] for x in HO_mut_exp if x[2]>0]
        cd_par_exp = [x[2] for x in CD_mut_exp if x[2]>0]
        
#         print('parallel')
#         print(sum(ho_par_exp))
#         print(len(ho_par_exp))
#         print(len(HO_mut_exp))
#         print()
#         print(sum(cd_par_exp))
#         print(len(cd_par_exp))
#         print(len(CD_mut_exp))
#         print()
#         print(sum(ho_par_exp)/sum(cd_par_exp))
#         quit()
        
        ratio_sites_par_exp = sum(ho_par_exp)/sum(cd_par_exp)
        ratio_genes_par_exp = len(ho_par_exp)/len(cd_par_exp)
        
        data_sim_par = np.genfromtxt('sum_parmut%s.csv' % suffix, dtype=None, 
                                     delimiter=',', names=True, encoding='UTF-8')
        
        sum_par_HO = [d[0] for d in data_sim_par]
        sum_par_CD = [d[1] for d in data_sim_par]
#         ratio_par_sim = [d[0]/d[1] for d in data_sim_par if d[1]>0]

#         plt.figure('R_N parallel')
#         dx = 0.1
#         bins = np.arange(0, 2.5+dx, dx)
#         plt.hist(ratio_par_sim, bins=bins, density=False, 
#                  weights=[1/len(ratio_par_sim)]*len(ratio_par_sim), 
#                  color='white', ec=ec, fc=fc, lw=lw, label=label)
#         plt.plot([ratio_sites_par_exp], [0.2], '*', ms=12, mec=fc, mfc=ec, label='observed')
#         plt.xlabel(r'$R_N$')
#         plt.ylabel('frequency')
# #         plt.xticks([bins[i] for i in range(len(bins)) if np.mod(i,4)==0])
#         plt.ylim(ymax=0.5)
#         plt.legend(loc=0)
        
        plt.figure('parmut%s' % suffix)
        bins = get_bins(sum_par_HO, sum_par_CD)
        plt.hist(sum_par_HO, bins=bins, color='white', fc='0.75', 
                 ec='red', lw=1, label='simulated HO')
        plt.hist(sum_par_CD, bins=bins, color='white', fc='None', 
                 ec='blue', lw=2, label='simulated CD')
        plt.annotate('Observed HO = {:,}'.format(sum(ho_par_exp)) , xy=(0.6,0.75), 
                     xycoords='axes fraction', weight='bold', color='red')
        plt.annotate('Observed CD = {:,}'.format(sum(cd_par_exp)), xy=(0.6,0.69), 
                     xycoords='axes fraction', weight='bold', color='blue')
        plt.xlabel('# parallel mutations') 
        plt.ylabel('number')
        # plt.xticks(range(0,11))
        plt.legend(loc=0)
    
    plt.show()
# '''
SUFFIX = ['_5strains', '_50strains']
EC = ['black', 'blue']
FC = ['0.75','None']
LW = [1,2]
LABEL = ['5 strains', '50 strains']

make_plots(SUFFIX, EC, FC, LW, LABEL)

quit()
# '''
########## RUN NEW SIMULATION ##########
filename = 'INPUT/gene_length_mut_5strains.csv'
suffix = re.match(".+mut(.*)\.csv",filename).group(1)
data = np.genfromtxt(filename, dtype=None, delimiter=',', names=True, encoding='UTF-8')

suffix += '_test'

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
    sum_hsmut_sites = [0,0]
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

with open('sum_hsmut_sites%s.csv' % suffix, 'w') as f:
    f.write("HO,CD\n")
    for n in range(niter):
        f.write("%d,%d\n" % (sum_hsmut_sites_HO[n],sum_hsmut_sites_CD[n]))

with open('sum_hsmut_genes%s.csv' % suffix, 'w') as f:
    f.write("HO,CD\n")
    for n in range(niter):
        f.write("%d,%d\n" % (sum_hsmut_genes_HO[n],sum_hsmut_genes_CD[n]))

with open('sum_parmut%s.csv' % suffix, 'w') as f:
    f.write("HO,CD\n")
    for n in range(niter):
        f.write("%d,%d\n" % (sum_parmut_HO[n],sum_parmut_CD[n]))

# plt.figure()
# plt.hist(ratio_sites, density=True, color='white', ec='black', lw=1)
# plt.xlabel(r'$R_N$')
# plt.ylabel('frequency')
# 
# plt.figure()
# plt.hist(ratio_genes, density=True, color='white', ec='black', lw=1)
# plt.xlabel(r'$R_G$')
# plt.ylabel('frequency')
# 
# plt.figure()
# plt.hist(sum_par_HO, density=False, color='white', fc='0.75', ec='red', lw=1, label='simulated HO')
# plt.hist(sum_par_CD, density=False, color='white', fc='None', ec='blue', lw=1, label='simulated CD')
# plt.annotate('Observed HO = %d' % 0  , xy=(0.6,0.75), xycoords='axes fraction', weight='bold', color='red')
# plt.annotate('Observed CD = %d' % 0 , xy=(0.6,0.69), xycoords='axes fraction', weight='bold', color='blue')
# plt.xlabel('# parallel mutations')
# plt.ylabel('frequency')
# plt.legend(loc=0)
# 
# plt.show()





