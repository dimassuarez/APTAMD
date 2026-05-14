#!/usr/bin/env python3
#

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Load the data
A = np.loadtxt('s1_plot.tabDUMMY_SUFFIX_TAB')
s1 = A[:, 1]
isnap = A[:, 0]
del A

cutoff = [ DUMMY_CUTOFF ]
nplot = len(cutoff)

B = np.loadtxt('DUMMY_METHOD_plot.tabDUMMY_SUFFIX_TAB')
n, m = B.shape

ncorr = m // nplot
j = 0
scc = np.zeros((n, nplot))
for i in range(ncorr, m + 1, ncorr):
    scc[:, j] = B[:, i - 1]
    j += 1

T = DUMMY_TEMP / 1000

# Plotting
plt.figure(1)
plt.plot(isnap, -T * s1, '-o', linewidth=2, markersize=4)
ymin = np.min(-T * s1)
ymax = np.max(-T * s1)
txt = ['S_1']
composite =  DUMMY_COMPOSITE

for j in range(nplot):
    plt.plot(isnap, -T * scc[:, j], '-o', linewidth=2, markersize=4)
    if composite == 1:
        txt.append('DUMMY_METHOD Comp. r_c= %5.2f'%(cutoff[j]))
    else:
        txt.append('DUMMY_METHOD r_c= %5.2f'%(cutoff[j]))
    ymin_corr = np.min(-T * scc[:, j])
    ymax_corr = np.max(-T * scc[:, j])
    ymin = min(ymin, ymin_corr)
    ymax = max(ymax, ymax_corr)

# Splines fitting of the limiting values
ns = 0
scc_lim = []
r_lim = []

for j in range(nplot):
    if cutoff[j] > -1:
        ns += 1
        scc_lim.append(scc[-1, j])
        r_lim.append(cutoff[j])

if ns > 0:
    r_lim = np.array(r_lim)
    scc_lim = np.array(scc_lim)
    rx = np.linspace(min(r_lim), max(r_lim), int((max(r_lim) - min(r_lim)) / 0.05) + 1)
    if ( len ( r_lim ) > 3 ):
        spkind='cubic'
    else:
        spkind='linear'
    spline = interp1d(r_lim, scc_lim, kind=spkind)
    scc_y = spline(rx)
    scc_opt = np.min(scc_y)
    imin = np.argmin(scc_y)

plt.grid(True)
plt.gcf().set_size_inches(8, 6)
plt.gca().set_facecolor('white')
plt.xlim([0, np.max(isnap)])
plt.ylim([ymin * 1.1, ymax * 0.9])
plt.xlabel('# Snap', fontsize=14)
plt.ylabel('-TS$_{conform}$ kcal/mol', fontsize=14)
plt.legend(txt, loc='best')
plt.savefig('s_DUMMY_METHOD_plotDUMMY_SUFFIX_TAB.png', dpi=300)

if ns > 0:
    plt.figure(2)
    plt.plot(rx, scc_y, '-r', linewidth=2)
    plt.plot(r_lim, scc_lim, 'ro', markersize=10)
    plt.plot(rx[imin], scc_opt, 'mo', markersize=12)
    ymin = np.min(scc_y)
    ymax = np.max(scc_y)
    plt.grid(True)
    plt.xlim([0.8 * min(r_lim), 1.1 * max(r_lim)])
    plt.ylim([ymin * 0.9, ymax * 1.1])
    plt.xlabel('Cutoff', fontsize=14)
    plt.ylabel('$S_{conform}^{limit}$ kcal/mol', fontsize=14)
    plt.title('S vs Cutoff:  $S_{min}$ = %6.2f  cal/mol K $-T S_{min}$= %6.2f kcal/mol at r= %5.2f Ang'%(scc_opt,-T*scc_opt,rx[imin]),fontsize=12)
    plt.legend(['Spline Fit', 'Calc S', 'Optimal'], loc='best')
    plt.gcf().set_size_inches(8, 6)
    plt.gca().set_facecolor('white')
    plt.savefig('s_DUMMY_METHOD_cutoffDUMMY_SUFFIX_TAB.png' if composite == 1 else 's_DUMMY_METHOD_cutoffDUMMY_SUFFIX_TAB.png', dpi=300)

# Write to file
filename = 's_DUMMY_METHOD_composite.tabDUMMY_SUFFIX_TAB' if composite == 1 else 's_DUMMY_METHOD.tabDUMMY_SUFFIX_TAB'
with open(filename, 'w') as fid:
    header = '# NumSnap , '
    for c in cutoff:
          if composite == 1 :
                header=header+ ' DUMMY_METHOD_Comp_rc= %5.2f ,'%c
          else:
                header=header+ ' DUMMY_METHOD_rc= %5.2f , '%c
    fid.write(header)
    fid.write('\n')
    for i in range(n):
        dataline=str( int( isnap[i] ) ) +', '
        for ic in  range(len(cutoff)):
              dataline=dataline+ ' %5.2f ,'%scc[i, ic]   
        fid.write(dataline)
        fid.write('\n')

