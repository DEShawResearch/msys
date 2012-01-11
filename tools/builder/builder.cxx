#include "builder.hxx"
#include "geom.hxx"
#include <boost/foreach.hpp>
#include <stdexcept>


namespace desres { namespace msys { namespace builder {

    namespace {
        Id find_atom(SystemPtr mol, Id res, Id prev, Id next, adef_t const& def) {
            if (fabs(def.rel)>1) {
                throw std::runtime_error(
                        "unsupported relative bond with magnitude > 1");
            }
            res = def.rel== 1 ? next :
                  def.rel==-1 ? prev :
                  res;
            if (!bad(res)) {
                BOOST_FOREACH(Id atm, mol->atomsForResidue(res)) {
                    if (def.name == mol->atom(atm).name) return atm;
                }
            }
            return BadId;
        }

        static const double def_bond = 1.0;
        static const double def_angle = 109.5 * M_PI/180;
        static const double def_dihedral = -120.0 * M_PI/180;

        bool wild_guess(SystemPtr mol, Id atm, IdSet const& added) {

            char name[256];
            sprintf(name, "%s %s %d", mol->atom(atm).name.c_str(), 
                    mol->residue(mol->atom(atm).residue).name.c_str(), 
                    mol->residue(mol->atom(atm).residue).resid);

            /* require a single bond to this atom */
            if (mol->bondCountForAtom(atm)!=1) {
                printf("Not building %s: bond count is %d\n",name,
                        mol->bondCountForAtom(atm));
                return false;
            }

            /* ensure the parent atom is sane */
            Id root = mol->bondedAtoms(atm)[0];
            IdList const& bonded = mol->bondedAtoms(root);
            if (bonded.size()>4) return false;

            /* find the identity of the other atoms bonded to root */
            Id a1=BadId;
            Id a2=BadId;
            Id a3=BadId;

            BOOST_FOREACH(Id b, bonded) {
                if (b==atm) continue;
                if (added.count(b)) continue;
                if      (bad(a2)) a2=b;
                else if (bad(a1)) a1=b;
                else if (bad(a3)) a3=b;
                else break;
            }

            Vec3 rtpos(          &mol->atom(root).x                   );
            Vec3 a1pos(bad(a1) ? &mol->atom(root).x : &mol->atom(a1).x);
            Vec3 a2pos(bad(a2) ? &mol->atom(root).x : &mol->atom(a2).x);
            Vec3 a3pos(bad(a3) ? &mol->atom(root).x : &mol->atom(a3).x);

            Float angle = def_angle;
            Float dihedral = def_dihedral;

            if (bonded.size()==3 && 
                    !bad(a1) && !bad(a2) && 
                    mol->atom(a1).atomic_number>1 &&
                    mol->atom(a2).atomic_number>1) { 

                angle = calc_angle(a1pos, rtpos, a2pos);
                angle = M_PI - 0.5*angle;
                dihedral = M_PI;

            } else if (bonded.size()==4 && 
                    !bad(a1) && !bad(a2) && !bad(a3)) {
                /* we have all but one of the atoms in a tetrahedral
                 * geometry.  find the position of the fourth by taking
                 * the plane formed by the non-root atoms */
                Vec3 v1(a1pos);
                v1 -= a2pos;
                Vec3 v2(a3pos);
                v2 -= a2pos;
                Vec3 d(rtpos);
                d -= a1pos;
                Vec3 v3 = calc_cross_prod(v1,v2);
                calc_normal(v3);
                v3 *= def_bond;

                if (d.dot(v3) <0) v3 *= -1;
                mol->atom(atm).x = rtpos.x + v3.x;
                mol->atom(atm).y = rtpos.y + v3.y;
                mol->atom(atm).z = rtpos.z + v3.z;
                return true;
            }

            Vec3 pos = apply_dihedral_geometry(a1pos, a2pos, rtpos,
                    def_bond, angle, dihedral);

            mol->atom(atm).x = pos.x;
            mol->atom(atm).y = pos.y;
            mol->atom(atm).z = pos.z;
            return true;
        }
    }

    void build( defs_t const& defs, SystemPtr mol, Id chain,
                std::string pfirst, std::string plast) {
    
        std::vector<resdef_t> resdefs;
        IdList const& residues = mol->residuesForChain(chain);
        IdSet added;
    
        printf("building chain '%s' with %d residues\n",
                mol->chain(chain).name.c_str(), (int)residues.size());

        /* assign resdef to each residue */
        for (Id i=0; i<residues.size(); i++) {
            Id ires=residues[i];
            residue_t& res = mol->residue(ires);
            ResdefMap::const_iterator idef = defs.resdefs.find(res.name);
            if (idef==defs.resdefs.end()) {
                fprintf(stderr, "No topology for residue %s:%d\n",
                        res.name.c_str(), res.resid);
                return;
            }
            resdef_t resdef = idef->second;

            /* remove all bonds in residue */
            BOOST_FOREACH(Id atm, mol->atomsForResidue(ires)) {
                mol->delBonds(mol->bondsForAtom(atm));
            }

            /* apply patch to resdef if first or last.  Hierarchy is
             * 1) arguments to function
             * 2) residue-specific patches
             * 3) global defaults in def_t.
             */
            if (i==0) {
                if (!pfirst.size()) pfirst = resdef.pfirst;
                if (!pfirst.size()) pfirst = defs.pfirst;
                if (pfirst.size() && strcasecmp(pfirst.c_str(), "none")) {
                    printf("Applying patch '%s' to %s:%d at start of chain '%s'\n",
                            pfirst.c_str(), 
                            res.name.c_str(),
                            res.resid,
                            mol->chain(chain).name.c_str());

                    ResdefMap::const_iterator pdef = defs.resdefs.find(pfirst);
                    if (pdef==defs.resdefs.end()) {
                        fprintf(stderr, "No such patch residue '%s'\n",
                                pfirst.c_str());
                        return;
                    }
                    pdef->second.patch_topology(resdef);
                }
            }
            if (i==residues.size()-1) {
                if (!plast.size()) plast = resdef.plast;
                if (!plast.size()) plast = defs.plast;
                if (plast.size() && strcasecmp(plast.c_str(), "none")) {
                    printf("Applying patch '%s' to %s:%d at end of chain '%s'\n",
                            plast.c_str(), 
                            res.name.c_str(),
                            res.resid,
                            mol->chain(chain).name.c_str());

                    ResdefMap::const_iterator pdef = defs.resdefs.find(plast);
                    if (pdef==defs.resdefs.end()) {
                        fprintf(stderr, "No such patch residue '%s'\n",
                                plast.c_str());
                        return;
                    }
                    pdef->second.patch_topology(resdef);
                }
            }

            /* Remove unrepresented atoms */
            IdList atoms = mol->atomsForResidue(ires);
            std::map<std::string, Id> amap;
            for (Id j=0; j<atoms.size(); j++) {
                msys::atom_t& atm = mol->atom(atoms[j]);
                AtomMap::const_iterator adef=resdef.atoms.find(atm.name);
                /* remove if not in resdef or if duplicate name */
                if (adef==resdef.atoms.end() || amap.count(atm.name)) {
                    //printf("deleted atom %s\n", atm.name.c_str());
                    mol->delAtom(atoms[j]);
                } else {
                    amap[atm.name] = atoms[j];
                }
            }
            /* Add missing atoms */
            BOOST_FOREACH(AtomMap::value_type const& adef, resdef.atoms) {
                Id atm = BadId;
                if (!amap.count(adef.first)) {
                    //printf("added atom %s\n", adef.first.c_str());
                    atm = mol->addAtom(ires);
                    added.insert(atm);
                } else {
                    atm = amap[adef.first];
                }
                mol->atom(atm).name = adef.first;
                mol->atom(atm).charge = adef.second.charge;

                TypeMap::const_iterator tdef = defs.types.find(adef.second.type);
                if (tdef==defs.types.end()) {
                    fprintf(stderr, "Invalid type '%s' for atom %s\n",
                            adef.second.type.c_str(), adef.first.c_str());
                } else {
                    mol->atom(atm).atomic_number = tdef->second.anum;
                    mol->atom(atm).mass = tdef->second.mass;
                }
            }

            /* cache resdef for later */
            resdefs.push_back(resdef);
        }

        /* add bonds.  We to have added all atoms before we can do this,
         * hence the need for a second pass. */
        for (Id i=0; i<residues.size(); i++) {
            resdef_t const& resdef = resdefs[i];
            Id res=residues[i];
            Id prev = i==0 ? BadId : residues[i-1];
            Id next = i==residues.size()-1 ? BadId : residues[i+1];

            BOOST_FOREACH(bond_t const& b, resdef.bonds) {
                Id atm1 = find_atom(mol, res, prev, next, b.def1);
                Id atm2 = find_atom(mol, res, prev, next, b.def2);
                if (bad(atm1) || bad(atm2)) {
                    printf("Skipping bond to nonexistent atom %s-%s\n",
                            b.def1.name.c_str(), b.def2.name.c_str());
                    continue;
                }
                mol->addBond(atm1, atm2);
            }
        }

        /* Add coordinates from conf records */
        for (Id i=0; i<residues.size(); i++) {
            resdef_t const& resdef = resdefs[i];
            Id res=residues[i];
            Id prev = i==0 ? BadId : residues[i-1];
            Id next = i==residues.size()-1 ? BadId : residues[i+1];

            BOOST_FOREACH(conf_t const& c, resdef.confs) {
                Id A = find_atom(mol, res, prev, next, c.def1);
                Id B = find_atom(mol, res, prev, next, c.def2);
                Id C = find_atom(mol, res, prev, next, c.def3);
                Id D = find_atom(mol, res, prev, next, c.def4);
                if (bad(A) || bad(B) || bad(C) || bad(D)) {
                    printf("Skipping geometry for %s:%d %s%s:%s:%s:%s%s\n",
                            mol->residue(res).name.c_str(),
                            mol->residue(res).resid,
                            c.def1.rel>0 ? "+" : c.def1.rel<0 ? "-" : "",
                            c.def1.name.c_str(), c.def2.name.c_str(),
                            c.def3.name.c_str(), 
                            c.def4.rel>0 ? "+" : c.def4.rel<0 ? "-" : "",
                            c.def4.name.c_str());
                    continue;
                }
                if (added.count(B) || added.count(C)) {
                    continue;
                }
                Float r=0, theta=0, phi=0;
                Vec3 apos, bpos, cpos;
                if ( mol->atom(D).atomic_number>1 && 
                     added.count(A)==0 && added.count(D)==1) { 
                    r=c.r2;
                    theta=c.a2;
                    phi=c.phi;
                    apos = Vec3(&mol->atom(A).x);
                    bpos = Vec3(&mol->atom(B).x);
                    cpos = Vec3(&mol->atom(C).x);
                } else if ( mol->atom(A).atomic_number>1 &&
                     added.count(A)==1 && added.count(D)==0) {
                    r=c.r1;
                    theta=c.a1;
                    phi=c.phi;
                    apos = Vec3(&mol->atom(D).x);
                    bpos = Vec3(&mol->atom(C).x);
                    cpos = Vec3(&mol->atom(B).x);
                    D=A;
                } else {
                    continue;
                }
                if (r==0) r=def_bond;
                if (theta==0) theta=def_angle;
                Vec3 dpos = apply_dihedral_geometry(apos,bpos,cpos,r,theta,phi);
                mol->atom(D).x=dpos.x;
                mol->atom(D).y=dpos.y;
                mol->atom(D).z=dpos.z;
                added.erase(added.find(D));

                /* if we guess a heavy atom position, guess its attached H */
                IdList bonded = mol->bondedAtoms(D);
                BOOST_FOREACH(Id b, bonded) {
                    if (mol->atom(b).atomic_number==1) {
                        added.insert(b);
                    }
                }
            } /* for each conf */
        } /* for each residue */

        /* make a chemically reasonable guess for atoms with one bond */
        IdSet tmp(added);
        BOOST_FOREACH(Id atm, tmp) {
            if (wild_guess(mol, atm, added)) {
                added.erase(added.find(atm));
            }
        }
        if (added.size()) printf("failed to guess positions for %d atoms\n",
                (int)added.size());

    }

    void patch( defs_t const& defs, std::string const& pres,
                SystemPtr mol, IdList const& residues ) {

        printf("applying patch %s\n", pres.c_str());
        ResdefMap::const_iterator it=defs.resdefs.find(pres);
        if (it==defs.resdefs.end()) {
            throw std::runtime_error("patch: no such patch definition");
        }
        resdef_t const& def = it->second;
        if (!def.patch) {
            throw std::runtime_error("patch: specified pres is not a patch");
        }
        if (residues.size()<1) {
            throw std::runtime_error("patch: need at least one residue");
        }

        /* remove atoms */
        for (unsigned i=0; i<def.delatoms.size(); i++) {
            adef_t const& adef = def.delatoms[i];
            if (adef.res >= (int)residues.size()) {
                fprintf(stderr, "patch residue %s requests residue %d, but only %d targets provided\n",
                        pres.c_str(), adef.res+1, (int)residues.size());
                throw std::runtime_error("patch: invalid targets");
            }
            assert(adef.res>=0);
            assert(adef.rel==0);    /* does this ever happen? */
            Id atm = find_atom(mol, residues.at(adef.res), BadId, BadId, adef);
            if (bad(atm)) {
                fprintf(stderr, "patch residue %s could find atom to delete: %s\n", pres.c_str(), adef.name.c_str());
                throw std::runtime_error("patch: bad target");
            }
            printf("deleting atom %s in residue %s:%d chain '%s'\n",
                    adef.name.c_str(), 
                    mol->residue(residues.at(adef.res)).name.c_str(),
                    mol->residue(residues.at(adef.res)).resid,
                    mol->chain(mol->residue(residues.at(adef.res)).chain).name.c_str());

            mol->delAtom(atm);
        }

        /* remove bonds? */
        assert (def.delbonds.size()==0);

        /* add/replace atoms */
        for (AtomMap::const_iterator it=def.atoms.begin(); it!=def.atoms.end();
                ++it) {

        }



    }
}}}
