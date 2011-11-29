

def CtBonds(e):
    ''' return a 0-based neighbor list for each atom. '''
    return [sorted(e.bondedAtoms(a)) for a in e.atoms()]

def stage2(molC, molB, atommaps, bmap,
           bondmaps, anglemaps, dihedmaps, 
           keep_extra_dihedral):

  # now we have the more complicated case of picking extra kept terms by 
  # choosing up to six restraints for blobs of dummies connected to blobs 
  # of real atoms.
  # First, construct a dictionary into anglemaps and dihedmaps so that we can
  # access individual terms and change their 0 to -1 (if kept)

  bond_dict = dict()
  for i,[[ti,tj],bond] in enumerate(bondmaps):
    bond_dict[tuple(bond)] = i
    bond_dict[tuple(reversed(bond))] = i
    
  angle_dict = dict()
  for i,[[ti,tj],angle] in enumerate(anglemaps):
    angle_dict[tuple(angle)] = i
    angle_dict[tuple(reversed(angle))] = i

  dihed_dict = dict()
  for i,[[ti,tj],dihed] in enumerate(dihedmaps):
    dihed_dict[tuple(dihed)] = i
    dihed_dict[tuple(reversed(dihed))] = i

  # first consider A state

  # contains 1-based atom numbers
  # backward edges are present
  neighbors = CtBonds(molC)

  # boolean list that indicates what atoms are always real and what 
  # atoms map to dummies in the B state
  realdummy = [0]*len(atommaps) # more than necessary
  for ai,aj in atommaps:
    if ai >= 0:
      if aj >= 0:
        realdummy[ai] = True
      else:
        realdummy[ai] = False

  # boolean list that indicates whether a dummy atom is connected through
  # kept covalent bonds to another dummy atom that is connected to a real
  # atom through a kept covalent bond.
  connected = [False]*len(atommaps)
  
  # loop through all bonds - only use (i,j) when i becomes a dummy in state B
  for i,jlist in enumerate(neighbors):
    if atommaps[i][0] < 0: continue
    if (not realdummy[i] and not connected[i]):
      for j in jlist:
        if realdummy[j]:
          # j is real and i is dummy in B state
          n1 = two_neighbors(i, neighbors, realdummy)
          n2 = two_neighbors(j, neighbors, realdummy)

          #print 'Stage A: Using template', n1, n2

          # Keep one bond that connects the dummy atom to the real atoms.
          ind = bond_dict.get( (n1[0],n2[0]), None )
          if (ind): keep(bondmaps, ind)
          
          # now we have five cases
          if len(n1) >= 2:
            ind = angle_dict.get( (n1[1],n1[0],n2[0]), None )
            if (ind): keep(anglemaps, ind)
          if len(n1) >= 3:
            ind = dihed_dict.get( (n1[2],n1[1],n1[0],n2[0]), None )
            if (ind and keep_extra_dihedral): keep(dihedmaps, ind)
          if len(n2) >= 2:
            ind = angle_dict.get( (n2[1],n2[0],n1[0]), None )
            if (ind): keep(anglemaps, ind)
          if len(n2) >= 3:
            ind = dihed_dict.get( (n2[2],n2[1],n2[0],n1[0]), None )
            if (ind and keep_extra_dihedral): keep(dihedmaps, ind)
          if len(n1) >= 2 and len(n2) >= 2:
            ind = dihed_dict.get( (n1[1],n1[0],n2[0],n2[1]), None )
            if (ind and keep_extra_dihedral): keep(dihedmaps, ind)
          break
      find_connected_dummies( i, neighbors, realdummy, connected )
      
  # now consider B state (refactor later)

  # contains 1-based atom numbers
  # backward edges are present
  neighbors = CtBonds(molB)

  # boolean list that indicates what atoms are always real and what 
  # atoms map to dummies in the A state
  realdummy = [0]*len(atommaps) # more than necessary
  for ai,aj in atommaps:
    if aj >= 0:
      if ai >= 0:
        realdummy[aj] = True
      else:
        realdummy[aj] = False

  # This bookkeeps whether a dummy atom is already connected - directly
  # or indirectly - to some real atom through a path of kept covalent
  # bond. Each group of connected dummy atoms should have one and only
  # one kept covalent bond to be connected to the real atoms.
  connected = [False]*len(atommaps)

  # loop through all bonds - only use (i,j) when i<j
  for i,jlist in enumerate(neighbors):
    if atommaps[i][0] < 0: continue
    # Seek a connection between the dummy atom and the real atoms, if and
    # only if it is not yet connected - directly or indirectly - to some
    # real atoms.
    if not realdummy[i] and not connected[i]:
      for j in jlist:
        if realdummy[j]:
          # j is real and i is dummy in A state
          n1 = two_neighbors(i, neighbors, realdummy)
          n2 = two_neighbors(j, neighbors, realdummy)

          # convert to combined numbering
          n1 = [bmap[x] for x in n1]
          n2 = [bmap[x] for x in n2]

          #print 'Stage B: Using template', n1, n2

          ind = bond_dict.get( (n1[0],n2[0]), None )
          if (ind): keep(bondmaps, ind)

          # now we have five cases
          if len(n1) >= 2:
            ind = angle_dict.get( (n1[1],n1[0],n2[0]), None )
            if (ind): keep(anglemaps, ind)
          if len(n1) >= 3:
            ind = dihed_dict.get( (n1[2],n1[1],n1[0],n2[0]), None )
            if (ind and keep_extra_dihedral): keep(dihedmaps, ind)
          if len(n2) >= 2:
            ind = angle_dict.get( (n2[1],n2[0],n1[0]), None )
            if (ind): keep(anglemaps, ind)
          if len(n2) >= 3:
            ind = dihed_dict.get( (n2[2],n2[1],n2[0],n1[0]), None )
            if (ind and keep_extra_dihedral): keep(dihedmaps, ind)
          if len(n1) >= 2 and len(n2) >= 2:
            ind = dihed_dict.get( (n1[1],n1[0],n2[0],n2[1]), None )
            if (ind and keep_extra_dihedral): keep(dihedmaps, ind)
          break
      find_connected_dummies( i, neighbors, realdummy, connected )

def find_connected_dummies( i, neighbors, realdummy, connected ):
  """
  Find all the dummy atoms that are connected to i - directly or
  indirectly, through covalent bonds between dummy atoms, and set
  the corresponding connected bits to true for these dummy atoms.
  """
  # We use a Depth-First-Search to find the connected component
  queue = [ i ]
  visited = [False]*len(realdummy)
  while len(queue) > 0:
    k = queue.pop()
    connected[k] = True
    visited[k] = True
    jlist = neighbors[k]
    for jp in jlist:
      j = jp 
      if (visited[j]): continue
      if (not realdummy[j]): continue
      queue.append( j )
  
def keep(termmaps, i):
  [ind,term] = termmaps[i]
  assert ind[0] == 0 or ind[1] == 0
  if ind[0] != 0 and ind[1] != 0:
    #print '*****Warning:', termmaps[i]
    return
  if ind[0] == 0:
    ind[0] = -1
  if ind[1] == 0:
    ind[1] = -1
  #print 'Keeping', term

def two_neighbors(i, neighbors, realdummy):
  '''Return a list of 1, 2, or 3 atoms [i, j, k] such that i is bonded to j, 
     j is bonded to k, and k is not i.  Input uses 0-based indexing.
     Neighbors list has both directions.
  '''
  ret = list()
  ret.append(i)

  # select neighbors that have the same type
  neigh = [ j for j in neighbors[i] if realdummy[j] == realdummy[i] ]
  if len(neigh) == 0:
    return ret  # could not find a j

  # here we know there is at least one j that works
  # pick the one with the largest degree (calculated using correct type)

  neigh_degree = list()
  for x in neigh:
    # y is the list of neighbors of x that have the right type
    y = [ z for z in neighbors[x] if realdummy[z] == realdummy[i] ]
    neigh_degree.append(len(y))

  # now choose the one with the largest degree
  ind = neigh_degree.index(max(neigh_degree))
  j = neigh[ind]
  ret.append(j)

  # now try to find second atom
  for k in neighbors[j]:
    if realdummy[k] == realdummy[i] and k != i:
      ret.append(k)
      break

  return ret


