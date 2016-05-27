#!/usr/bin/env python2.7

import sys

sys.path.insert(0,'objs/CentOS6/x86_64/lib/python')
import msys
import vmd
import molecule
import atomsel
from time import time

def compare_atomsel(coord_ent,sel, dump=False, perf=False):
  try:
    ent_gids = coord_ent.selectIds(sel)
  except RuntimeError:
    print "MSYS failed to parse %s" % sel
    return
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
    rc=False
  else:
    rc=True
    print "match for [%s] (%d) "% (sel, len(vmd_gids))
  if dump:
      print "MSYS:" 
      print ent_gids
      print "VMD:"
      print vmd_gids
  if perf:
      t0=time()
      coord_ent.selectIds(sel)
      t1=time()
      atomsel.atomsel(sel)
      t2=time()
      print "MSYS: %s ms" % ((t1-t0)*1000)
      print "VMD:  %s ms" % ((t2-t1)*1000)
  return rc

if len(sys.argv)<2:
    DMSROOT='/proj/desres/root/Linux/x86_64/dms_inputs/1.5.5/share'
    path='%s/2f4k.dms' % DMSROOT
else:
    path=sys.argv[1]

coord_vmd=molecule.load(path.split('.')[-1], path)
coord_ent = msys.Load(path)
coord_ent.select('none')

compare_atomsel(coord_ent,"protein or water and element O")
compare_atomsel(coord_ent,"protein and water or element O")
compare_atomsel(coord_ent,"protein or (water and element O)")
compare_atomsel(coord_ent,"(protein and water) or element O")

compare_atomsel(coord_ent,"ctnumber 1")
compare_atomsel(coord_ent,"ctnumber 2")
compare_atomsel(coord_ent,"ctnumber 3")
compare_atomsel(coord_ent,"ctnumber 4")
compare_atomsel(coord_ent,"all")
compare_atomsel(coord_ent,"none")
compare_atomsel(coord_ent,"index 10 20 30",perf=True)
compare_atomsel(coord_ent,"index 10 20 30",perf=True)
compare_atomsel(coord_ent,"index 10 20 30",perf=True)
compare_atomsel(coord_ent,"resid 1 to 20 22 24 to 27 35",perf=True)
compare_atomsel(coord_ent,"protein",perf=True)
compare_atomsel(coord_ent,"backbone",perf=True)
compare_atomsel(coord_ent,"sidechain",perf=True)
#compare_atomsel(coord_ent,"hydrogen")
compare_atomsel(coord_ent,"atomicnumber 1",perf=True)

compare_atomsel(coord_ent,"chain A B",perf=True)
#compare_atomsel(coord_ent,"segid L14")

compare_atomsel(coord_ent,"water",perf=True)
compare_atomsel(coord_ent,"x<3")
compare_atomsel(coord_ent,"y>=1")
compare_atomsel(coord_ent,"z<=-5")

compare_atomsel(coord_ent,"pbwithin 3 of protein", perf=True)
compare_atomsel(coord_ent,"(not atomicnumber 1) and pbwithin 3 of protein")
compare_atomsel(coord_ent,"backbone and pbwithin 4 of water", perf=True)

compare_atomsel(coord_ent,"same residue as within 4 of protein", perf=True)
compare_atomsel(coord_ent,"(not atomicnumber 1) and within 1.5 of protein", perf=True)
compare_atomsel(coord_ent,"(not atomicnumber 1) and within 2.5 of protein", perf=True)
compare_atomsel(coord_ent,"(not atomicnumber 1) and within 5.5 of protein", perf=True)
compare_atomsel(coord_ent,"within 0 of protein")
compare_atomsel(coord_ent,"withinbonds 2 of name C CA N O")
compare_atomsel(coord_ent,"fragment 0")
compare_atomsel(coord_ent,"fragment 1")
compare_atomsel(coord_ent,"not oxygen or water")
compare_atomsel(coord_ent,"(not oxygen) or water")
compare_atomsel(coord_ent,"not (oxygen or water)")
compare_atomsel(coord_ent,"not oxygen and water")
compare_atomsel(coord_ent,"(not oxygen) and water")
compare_atomsel(coord_ent,"not (oxygen and water)")

compare_atomsel(coord_ent,"sqr(x)/36 + sqr(y)/81 + sqr(z)/125 < 1")


compare_atomsel(coord_ent,"sqrt(sqr(x)+sqr(y))<5")
compare_atomsel(coord_ent,"abs(x-y)<5")

compare_atomsel(coord_ent,"residue % 10 == 0")
compare_atomsel(coord_ent,"residue^3 == 125")

# nearest
compare_atomsel(coord_ent,"water and noh")
compare_atomsel(coord_ent,"water and noh and nearest 1 to protein", perf=True)
compare_atomsel(coord_ent,"water and noh and nearest 10 to protein", perf=True)
compare_atomsel(coord_ent,"water and noh and nearest 100 to protein", perf=True)

compare_atomsel(coord_ent,"nearest 1 to protein", perf=True)
compare_atomsel(coord_ent,"nearest 10 to protein", perf=True)
compare_atomsel(coord_ent,"nearest 100 to protein", perf=True)
compare_atomsel(coord_ent,"nearest 1000 to protein", perf=True)
compare_atomsel(coord_ent,"nearest 20000 to protein", perf=True)

compare_atomsel(coord_ent,"water and same residue as nearest 1 to protein", perf=True)
compare_atomsel(coord_ent,"water and same residue as nearest 10 to protein", perf=True)
compare_atomsel(coord_ent,"water and same residue as nearest 100 to protein", perf=True)

# operator precedence

compare_atomsel(coord_ent,"x + y * z < 3")
compare_atomsel(coord_ent,"x * y - z < 3")


# test regex
compare_atomsel(coord_ent,'name "C.*"')
compare_atomsel(coord_ent,'name "C[a-z]"')
compare_atomsel(coord_ent,'name "C[A-Z]"')


# element
compare_atomsel(coord_ent,'element H')
compare_atomsel(coord_ent,'element C')
compare_atomsel(coord_ent,'element N')
compare_atomsel(coord_ent,'element O')
compare_atomsel(coord_ent,'element P')

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

