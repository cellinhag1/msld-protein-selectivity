import numpy as np
import os

alf_info = {}
alf_info['name']       = 'apo'
alf_info['nsubs']      = [2, 2, 2, 2, 2, 2, 2, 2]  # 8 INDEPENDENT sites x 2 states each
alf_info['nblocks']    = int(np.sum(alf_info['nsubs']))  # = 16
alf_info['ncentral']   = 0
alf_info['nreps']      = 1
alf_info['nnodes']     = 1
alf_info['enginepath'] = os.environ.get('CHARMM_EXE', 'charmm')  # path to the CHARMM binary
alf_info['temp']       = 298.15
alf_info['impcons']    = 'custom'   # pre-computed theta_bias_ind8.5 (fnex=8.5) via G_imp/ symlink
# Same 8 residues as apo_patch, now independent (own lambda each) instead of
# coupled into 1 shared lambda. Order matches alchemical_definitions.inp:
# 663(S->L) 778(I->V) 787(F->W) 804(L->M) 813(I->L) 816(M->L) 823(A->F) 824(I->V)
# All 8 swaps are charge-neutral: no DSC correction needed on any site.
alf_info['q'] = np.zeros(16, dtype=float)
