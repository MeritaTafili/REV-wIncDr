# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 19:21:41 2023

@author: merit
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
          
# -----------------------------------------------------------------------------
# Functions / Definitions
# -----------------------------------------------------------------------------
def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

file1 = 'CU_MCC_OCR1.out'
file2 = 'CU_MCC_OCR2.out'
file3 = 'CU_MCC_OCR4.out'
# -----------------------------------------------------------------------------
# Global settings
# -----------------------------------------------------------------------------
height = 6 # cm
width =  8.5  # cm
#%%
fig, (sp1,sp2) = plt.subplots(nrows=1, ncols=2, figsize=cm2inch(2*width,height), dpi=200, facecolor='white')

#
# PLOT 1
#

data = pd.read_csv(file1, skipinitialspace=False, delim_whitespace=True )
eps11 = -data['stran(1)']*100. ; eps22 = -data['stran(2)']*100. ; eps33 = -data['stran(3)']*100.
s11 = -(data['stress(1)']) ; s22 = -(data['stress(2)']) ; s33 = -(data['stress(3)'])
p = ( s11 + s22 + s33 ) / 3.
q = (s11-s22)

sp1.plot(eps11, q, color = 'blue', ls='-', lw=0.75, zorder=1, label='OCR=1')

data = pd.read_csv(file2, skipinitialspace=False, delim_whitespace=True )
eps11 = -data['stran(1)']*100. ; eps22 = -data['stran(2)']*100. ; eps33 = -data['stran(3)']*100.
s11 = -(data['stress(1)']) ; s22 = -(data['stress(2)']) ; s33 = -(data['stress(3)'])
p = ( s11 + s22 + s33 ) / 3.
q = (s11-s22)

sp1.plot(eps11, q, color = 'red', ls='-', lw=0.75, zorder=1, label='OCR=2')

data = pd.read_csv(file3, skipinitialspace=False, delim_whitespace=True )
eps11 = -data['stran(1)']*100. ; eps22 = -data['stran(2)']*100. ; eps33 = -data['stran(3)']*100.
s11 = -(data['stress(1)']) ; s22 = -(data['stress(2)']) ; s33 = -(data['stress(3)'])
p = ( s11 + s22 + s33 ) / 3.
q = (s11-s22)

sp1.plot(eps11, q, color = 'green', ls='-', lw=0.75, zorder=1, label='OCR=4')

# data = pd.read_csv(file4, skipinitialspace=False, delim_whitespace=True )
# eps11 = -data['stran(1)']*100. ; eps22 = -data['stran(2)']*100. ; eps33 = -data['stran(3)']*100.
# s11 = -(data['stress(1)']) ; s22 = -(data['stress(2)']) ; s33 = -(data['stress(3)'])
# p = ( s11 + s22 + s33 ) / 3.
# q = (s11-s22)

# sp1.plot(eps11, q, color = 'grey', ls='-', lw=0.75, zorder=1, label='he_s')

# sp1.set_xlim(-15.0,15.0)
# sp1.set_xticks(np.arange(-15,20,5))
sp1.set_ylim(0,450)
# sp1.set_yticks(np.arange(-500,600,100))
sp1.set_xlabel('$\\varepsilon_{1}$ in %')
sp1.set_ylabel('$q$ in kPa' )

# sp1.legend(loc='best')

ax = fig.gca()

#
# PLOT 2
#

p = np.arange(0.,550.,1)
Mc = 1.0*p
Me = 0.75*Mc
sp2.plot(p, Mc, 'r', lw=0.5, ls=':')
# sp2.plot(p, -Me, 'k', lw=0.5, ls=':')



data = pd.read_csv(file1, skipinitialspace=False, delim_whitespace=True )
eps11 = -data['stran(1)']*100. ; eps22 = -data['stran(2)']*100. ; eps33 = -data['stran(3)']*100.
s11 = -(data['stress(1)']) ; s22 = -(data['stress(2)']) ; s33 = -(data['stress(3)'])
p = ( s11 + s22 + s33 ) / 3.
q = (s11-s22)

sp2.plot(p, q, color = 'blue', ls='-', lw=0.75, zorder=1, label='OCR=1')

data = pd.read_csv(file2, skipinitialspace=False, delim_whitespace=True )
eps11 = -data['stran(1)']*100. ; eps22 = -data['stran(2)']*100. ; eps33 = -data['stran(3)']*100.
s11 = -(data['stress(1)']) ; s22 = -(data['stress(2)']) ; s33 = -(data['stress(3)'])
p = ( s11 + s22 + s33 ) / 3.
q = (s11-s22)

sp2.plot(p, q, color = 'red', ls='-', lw=0.75, zorder=1, label='OCR=2')
        
data = pd.read_csv(file3, skipinitialspace=False, delim_whitespace=True )
eps11 = -data['stran(1)']*100. ; eps22 = -data['stran(2)']*100. ; eps33 = -data['stran(3)']*100.
s11 = -(data['stress(1)']) ; s22 = -(data['stress(2)']) ; s33 = -(data['stress(3)'])
p = ( s11 + s22 + s33 ) / 3.
q = (s11-s22)

sp2.plot(p, q, color = 'green', ls='-', lw=0.75, zorder=1, label='OCR=4')  

# data = pd.read_csv(file4, skipinitialspace=False, delim_whitespace=True )
# eps11 = -data['stran(1)']*100. ; eps22 = -data['stran(2)']*100. ; eps33 = -data['stran(3)']*100.
# s11 = -(data['stress(1)']) ; s22 = -(data['stress(2)']) ; s33 = -(data['stress(3)'])
# p = ( s11 + s22 + s33 ) / 3.
# q = (s11-s22)

# sp2.plot(p, q, color = 'grey', ls='-', lw=0.75, zorder=1)  

sp2.set_xlim(0,450)
# sp2.set_xticks(np.arange(0,550,100))
sp2.set_ylim(0,450)
# sp2.set_yticks(np.arange(-500,600,100))
sp2.set_xlabel('$p$ in kPa')
sp2.set_ylabel('$q$ in kPa' )
sp2.legend(loc='best')
ax = fig.gca()


plt.tight_layout(w_pad=1.2)
plt.show()
fig.savefig('triax_CU.pdf', bbox_inches='tight')
fig.savefig('triax_CU.png', bbox_inches='tight')