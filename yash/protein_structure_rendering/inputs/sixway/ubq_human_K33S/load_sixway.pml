# 6-way comparison: ubq_human K33S (1UBQ vs 4PIH)
# exp = crystal (grey/black), af = AlphaFold2 (blue), esm = ESMFold (orange). WT light, mutant dark.
# Superposition: every object onto exp_wt, C-alpha only, cycles=0 (no outlier rejection,
# so real motion is not hidden). RMSDs print to the PyMOL console.
# Nothing here ray-traces or writes images.
# Run from inside this folder:  cd <this folder> && pymol load_sixway.pml
reinitialize
load exp_wt.pdb
load exp_mut.pdb
load af_wt.pdb
load af_mut.pdb
load esm_wt.pdb
load esm_mut.pdb
remove hydro
align exp_mut and name CA, exp_wt and name CA, cycles=0, object=aln_exp_mut
align af_wt and name CA, exp_wt and name CA, cycles=0, object=aln_af_wt
align af_mut and name CA, exp_wt and name CA, cycles=0, object=aln_af_mut
align esm_wt and name CA, exp_wt and name CA, cycles=0, object=aln_esm_wt
align esm_mut and name CA, exp_wt and name CA, cycles=0, object=aln_esm_mut
delete aln_*
hide everything
show cartoon
set cartoon_transparency, 0
color grey70, exp_wt
color grey30, exp_mut
color lightblue, af_wt
color tv_blue, af_mut
color lightorange, esm_wt
color orange, esm_mut
select site, resi 33
show sticks, site and not name N+C+O
util.cnc site
set grid_mode, 1
set grid_slot, 1, exp_wt
set grid_slot, 1, exp_mut
set grid_slot, 2, af_wt
set grid_slot, 2, af_mut
set grid_slot, 3, esm_wt
set grid_slot, 3, esm_mut
set grid_mode, 0
orient exp_wt
deselect
