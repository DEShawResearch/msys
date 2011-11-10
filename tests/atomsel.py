#!/usr/bin/env python2.7

import sys

sys.path.insert(0,'objs/Linux/x86_64/lib/python')
import msys
import vmd
import molecule
import atomsel

def compare_atomsel(coord_ent,sel):
  al = coord_ent.select(sel)
  ent_gids = []
  for a in al:
    ent_gids.append(a.id)

  vmd_atomsel = atomsel.atomsel(sel)
  vmd_gids = vmd_atomsel.get("index")

  if not ent_gids==vmd_gids:
    print "mismatch for [%s]: vmd %d msys %d" %(sel,len(vmd_gids),len(ent_gids))
    #print list(set(vmd_gids)-set(ent_gids))
    #print list(set(ent_gids)-set(vmd_gids))
    #print ent_gids
    #all=atomsel.atomsel()
    #names=all.get('name')
    #for g in ent_gids: 
      #if g not in vmd_gids:
        #print g,names[g], al[g].pos[:2]
    return False
  else:
    print "match for [%s] (%d) "% (sel, len(vmd_gids))
    return True

DMSROOT='/proj/desres/root/Linux/x86_64/dms_inputs/1.5.5/share'
coord_vmd=molecule.load("dms","%s/2f4k.dms" % DMSROOT)
coord_ent = msys.LoadDMS("%s/2f4k.dms" % DMSROOT,True)
#coord_vmd=molecule.load("dms","/tmp/pro.dms")
#coord_ent = msys.LoadDMS("/tmp/pro.dms",True)

compare_atomsel(coord_ent,"all")
compare_atomsel(coord_ent,"none")
compare_atomsel(coord_ent,"index 10 20 30")
compare_atomsel(coord_ent,"gid 10 20 30")
compare_atomsel(coord_ent,"resid 1 to 20 22 24 to 27 35")
compare_atomsel(coord_ent,"protein")
compare_atomsel(coord_ent,"backbone")
#compare_atomsel(coord_ent,"hydrogen")
compare_atomsel(coord_ent,"atomicnumber 1")

compare_atomsel(coord_ent,"chain A B")
#compare_atomsel(coord_ent,"segid L14")

compare_atomsel(coord_ent,"water")
compare_atomsel(coord_ent,"x<3")
compare_atomsel(coord_ent,"y>=1")
compare_atomsel(coord_ent,"z<=-5")

compare_atomsel(coord_ent,"pbwithin 3 of protein")
compare_atomsel(coord_ent,"(not atomicnumber 1) and pbwithin 3 of protein")
compare_atomsel(coord_ent,"backbone and pbwithin 4 of water")

compare_atomsel(coord_ent,"same residue as within 4 of protein")
compare_atomsel(coord_ent,"(not atomicnumber 1) and within 1.5 of protein")
compare_atomsel(coord_ent,"(not atomicnumber 1) and within 2.5 of protein")
compare_atomsel(coord_ent,"(not atomicnumber 1) and within 5.5 of protein")
compare_atomsel(coord_ent,"within 0 of protein")
compare_atomsel(coord_ent,"withinbonds 2 of name C CA N O")
compare_atomsel(coord_ent,"fragment 0")
compare_atomsel(coord_ent,"fragment 1")


compare_atomsel(coord_ent,"sqrt(sqr(x)+sqr(y))<5")
compare_atomsel(coord_ent,"abs(x-y)<5")

compare_atomsel(coord_ent,"residue % 10 == 0")


# test regex
compare_atomsel(coord_ent,'name "C.*"')
compare_atomsel(coord_ent,'name "C[a-z]"')
compare_atomsel(coord_ent,'name "C[A-Z]"')

# macros
for x in '''
at
acidic
cyclic
acyclic
aliphatic
alpha
amino
aromatic
basic
bonded
buried
cg
charged
hetero
hydrophobic
small
medium
large
neutral
polar
purine
pyrimidine
surface
lipid
lipids
ion
ions
sugar
solvent
carbon
hydrogen
nitrogen
oxygen
sulfur
noh
heme
'''.split():
    if x: 
        compare_atomsel(coord_ent,x)

