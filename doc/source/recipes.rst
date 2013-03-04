

Adding energy groups
--------------------

Desmond and Anton use the "energy_groups" atom property to assign atoms to
energy groups::

  mol = msys.Load('system.dms')

  # add an integer property.  The default value is zero.  It's a no-op
  # if the property already exists, and an error if it exists but has a
  # different type.
  mol.addAtomProp('grp_energy', int)        

  # assign protein to energy group 1
  for a in mol.select('protein'):
    a['grp_energy'] = 1

  # save the result
  msys.SaveDMS(mol, 'system_engrp.dms')



