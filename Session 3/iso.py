# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 19:21:41 2023

@author: merit
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from numpy import log as ln

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

file1 = 'output_iso_1_le.out'
file2 = 'output_iso_1_he_c.out'
file3 = 'output_iso_1_he_s.out'
file4 = 'output_iso_r_1_he_s.out'

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

sp1.semilogx(p, eps11, color = 'blue', ls='-', lw=0.75, zorder=1, label='le')

data = pd.read_csv(file2, skipinitialspace=False, delim_whitespace=True )
eps11 = -data['stran(1)']*100. ; eps22 = -data['stran(2)']*100. ; eps33 = -data['stran(3)']*100.
s11 = -(data['stress(1)']) ; s22 = -(data['stress(2)']) ; s33 = -(data['stress(3)'])
p = ( s11 + s22 + s33 ) / 3.
q = (s11-s22)

sp1.semilogx(p, eps11, color = 'red', ls='-', lw=0.75, zorder=1, label='he_c')

data = pd.read_csv(file3, skipinitialspace=False, delim_whitespace=True )
eps11 = -data['stran(1)']*100. ; eps22 = -data['stran(2)']*100. ; eps33 = -data['stran(3)']*100.
s11 = -(data['stress(1)']) ; s22 = -(data['stress(2)']) ; s33 = -(data['stress(3)'])
p = ( s11 + s22 + s33 ) / 3.
q = (s11-s22)

sp1.semilogx(p, eps11, color = 'green', ls='-', lw=0.75, zorder=1, label='he_s')

data = pd.read_csv(file4, skipinitialspace=False, delim_whitespace=True )
eps11 = -data['stran(1)']*100. ; eps22 = -data['stran(2)']*100. ; eps33 = -data['stran(3)']*100.
s11 = -(data['stress(1)']) ; s22 = -(data['stress(2)']) ; s33 = -(data['stress(3)'])
p = ( s11 + s22 + s33 ) / 3.
q = (s11-s22)

sp1.semilogx(p, eps11, color = 'grey', ls='-', lw=0.75, zorder=1, label='he_s_r')

# sp1.set_xlim(-15.0,15.0)
# sp1.set_xticks(np.arange(-15,20,5))
sp1.set_ylim(16,0)
# sp1.set_yticks(np.arange(-500,600,100))
sp1.set_ylabel('$\\varepsilon_{1}$ in %')
sp1.set_xlabel('ln$(p)$ in kPa' )

ax = fig.gca()

sp1.legend(loc='best')

#
# PLOT 2
#

p = np.arange(0.,550.,1)
Mc = 1.25*p
Me = 0.75*Mc
# sp2.plot(p, Mc, 'k', lw=0.5, ls=':')
# sp2.plot(p, -Me, 'k', lw=0.5, ls=':')



data = pd.read_csv(file1, skipinitialspace=False, delim_whitespace=True )
eps11 = -data['stran(1)']*100. ; eps22 = -data['stran(2)']*100. ; eps33 = -data['stran(3)']*100.
epsv = eps11 + eps22 + eps33
s11 = -(data['stress(1)']) ; s22 = -(data['stress(2)']) ; s33 = -(data['stress(3)'])  
epor = data['statev(']
p = ( s11 + s22 + s33 ) / 3.
q = (s11-s22)

# sp2.semilogx(p, epor, color = 'blue', ls='-', lw=0.75, zorder=1)

data = pd.read_csv(file2, skipinitialspace=False, delim_whitespace=True )
eps11 = -data['stran(1)']*100. ; eps22 = -data['stran(2)']*100. ; eps33 = -data['stran(3)']*100.
epsv = eps11 + eps22 + eps33
s11 = -(data['stress(1)']) ; s22 = -(data['stress(2)']) ; s33 = -(data['stress(3)'])
epor = data['statev(']
p = ( s11 + s22 + s33 ) / 3.
q = (s11-s22)

sp2.semilogx(p, epor, color = 'red', ls='-', lw=0.75, zorder=1, label='he_c')
kappa = 0.03
sp2.semilogx(p, 1.05-kappa*ln(p), color = 'red', ls='--', label='$e=e_0-\kappa*\ln(p)$', lw=0.75, zorder=1)

   
data = pd.read_csv(file3, skipinitialspace=False, delim_whitespace=True )
eps11 = -data['stran(1)']*100. ; eps22 = -data['stran(2)']*100. ; eps33 = -data['stran(3)']*100.
epsv = eps11 + eps22 + eps33
s11 = -(data['stress(1)']) ; s22 = -(data['stress(2)']) ; s33 = -(data['stress(3)'])
epor = data['statev(']
p = ( s11 + s22 + s33 ) / 3.
q = (s11-s22)

sp2.semilogx(p, epor, color = 'green', ls='-', lw=0.75, zorder=1, label='he_s')

data = pd.read_csv(file4, skipinitialspace=False, delim_whitespace=True )
eps11 = -data['stran(1)']*100. ; eps22 = -data['stran(2)']*100. ; eps33 = -data['stran(3)']*100.
epsv = eps11 + eps22 + eps33
s11 = -(data['stress(1)']) ; s22 = -(data['stress(2)']) ; s33 = -(data['stress(3)'])
epor = data['statev(']
p = ( s11 + s22 + s33 ) / 3.
q = (s11-s22)

sp2.semilogx(p, epor, color = 'grey', ls='-', lw=0.75, zorder=1, label='he_s_r')

# sp2.set_xlim(0,550)
# sp2.set_xticks(np.arange(0,550,100))
# sp2.set_ylim(-550,550)
# sp2.set_yticks(np.arange(-500,600,100))
sp2.set_xlabel('ln$(p)$ in kPa')
sp2.set_ylabel('$e$ ' )
sp2.legend(loc='best',handlelength=0.7)

ax = fig.gca()


plt.tight_layout(w_pad=1.2)
plt.show()
fig.savefig('iso.pdf', bbox_inches='tight')
fig.savefig('iso.png', bbox_inches='tight')