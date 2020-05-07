// @COPYRIGHT@
#include "prep_alchemical_mae.hxx"
#include "../../types.hxx"

#include <sstream>

#include "destro/Destro.hxx"

#include <algorithm>
#include <set>
#include <cctype>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <cstring>

namespace {
  /*! Provide case-insensitive string comparisons */
  struct CaseInsensitive {
    static std::string to_upper(const std::string& str) {
      std::string r(str);
      for (unsigned i=0; i<str.size(); i++) r[i] = std::toupper(str[i]);
      return r;
    }
    // for equality tests
    static bool equals( const std::string &s1, const std::string &s2 ) {
      std::string::const_iterator i1=s1.begin(), i2=s2.begin();
      for (; i1 != s1.end() && i2 != s2.end(); ++i1, ++i2) 
        if (std::toupper(*i1) != std::toupper(*i2)) return false;
      return s1.size() == s2.size();
    }
    // for associative containers
    bool operator()(const std::string &s1, const std::string &s2) const {
      std::string::const_iterator i1=s1.begin(), i2=s2.begin();
      for (; i1 != s1.end() && i2 != s2.end(); ++i1, ++i2) {
        char c1 = std::toupper(*i1);
        char c2 = std::toupper(*i2);
        if (c1<c2) return true;
        if (c2<c1) return false; 
      }
      return s1.size() < s2.size();
    }
  };
}

namespace {
  /*! entries in fepio_fep table */
  struct fep_elem {
    //! variable names correspond to fepio_fep field names
    //!@{
    int ti, tj;
    int ai, aj, ak, al;
    int am, an, ao, ap; // Allow treatment of torsion-torsion terms
    int moiety; //!< The moiety that each mapped term belongs to.
    //!@}
  };
  /*! Represnts a single table in fepio_fep */
  typedef std::vector<fep_elem> FepList;
  /*! Represnts the entire fepio_fep section */
  typedef std::map<std::string, FepList > FepioMapping;

  /*! create an FepList
   * @param[out] v fep table
   * @param[in]  m subsection of fepio_fep 
   */
  void parse_fepio_map( FepList &v, const desres::msys::Destro &m ) {
    for (unsigned i=0; i<m.size(); i++) {
      const desres::msys::Destro &record = m[i+1];
      fep_elem elem;
      elem.ti = record("fepio_ti").or_else(0);
      elem.tj = record("fepio_tj").or_else(0);
      // Make refer to combined ct atom numbers 
      elem.ai = record("fepio_ai").or_else(0);
      elem.aj = record("fepio_aj").or_else(0);
      elem.ak = record("fepio_ak").or_else(0);
      elem.al = record("fepio_al").or_else(0);
      elem.am = record("fepio_am").or_else(0);
      elem.an = record("fepio_an").or_else(0);
      elem.ao = record("fepio_ao").or_else(0);
      elem.ap = record("fepio_ap").or_else(0);
      // Get the moiety index
      elem.moiety = record("fepio_moiety").or_else(0);
      v.push_back(elem);
    }
  }

  /*! parse the entire fepio_fep block
   * @param[out] m parsed fepio_fep block
   * @param[in] fepio_fep section
   */
  void grovel_fepio_mapping( FepioMapping &m, const desres::msys::Destro &fepio ) {
    static const char *sections[] = {
      "fepio_atommaps", "fepio_bondmaps", "fepio_anglemaps",
      "fepio_dihedmaps", "fepio_pairmaps", "fepio_exclmaps",
      "fepio_morsebondmaps"
    };
    unsigned i,n = sizeof(sections)/sizeof(const char *);
    for (i=0; i<n; i++) {
      const char *s = sections[i];
      if (fepio.contains(s)) {
        FepList v;
        parse_fepio_map( v, fepio[s] );
        m[s] = v;
      }
    }
  }

  /*! Add DMY atom to vdw table if not already present
   * @param m1 ffio_ff section
   */
  void provide_dummy( desres::msys::Destro& m1) {
    desres::msys::Destro& vm1 = m1.block("ffio_vdwtypes");
    for ( unsigned int i = 1; i <= vm1->size(); ++i ) {
      if ( vm1[i]["ffio_name"] == "DMY" ) return;
    }
    desres::msys::Destro& nvm = vm1.append();
    nvm["ffio_name"] = "DMY";
    nvm["ffio_funct"] = "LJ12_6_sig_epsilon";
    nvm["ffio_c1"] = 0.0;   // sigma
    nvm["ffio_c2"] = 0.0;   // epsilon
  }

  /*! merge two force field entries into one alchemical entry
   * @param m1 alchemical stage 1 ffio_ff block
   * @param m2 alchemical stage 2 ffio_ff block
   * @param map iterator over fepio mapping
   * @param keeper field names that don't go to zero if target is zero
   * @param blk ffio_ff subsection
   */
  void combine_entry(desres::msys::Destro& m1, desres::msys::Destro const& m2,
      FepioMapping::const_iterator const& map, std::string const& keeper, 
      std::string const& blk, std::string const& alc_blk) {

    const char * moiety = "ffio_moiety";

    std::set<int> remove; // rows of m1[blk] to remove
    for (FepList::const_iterator i = map->second.begin(); i != map->second.end(); ++i) {
      char const* Cs[] = {
	"ffio_c0", "ffio_c1", "ffio_c2", "ffio_c3",
	"ffio_c4", "ffio_c5", "ffio_c6", "ffio_c7"
      };
      int const nCs = sizeof(Cs)/sizeof(char const*);
      int ta = i->ti;
      int tb = i->tj;
      if ( ta > 0 && tb == -1 ) {
        // Keep the A state intact, no alchemical transformation for the term.
        continue;

      } else if ( ta > 0 && tb == -2 ) {
        desres::msys::Destro &b1 = m1.block(alc_blk).append(m1[blk][ta]);
        b1.setattr( moiety, i->moiety );
        remove.insert(ta);
	for ( int j = 0; j < nCs; ++j ) {
          const std::string &Cj  = Cs[j];
          const std::string &CjB = Cj+'B';
          if (b1.has_value(Cj)) {
            b1.setattr(CjB, b1[Cj]);
          }
        }

      } else if ( ta > 0 && tb == 0 ) {
        // B state has zeroed coefficients
        desres::msys::Destro &b1 = m1.block(alc_blk).append(m1[blk][ta]);
        b1.setattr( moiety, i->moiety );
        remove.insert(ta);
	for ( int j = 0; j < nCs; ++j ) {
          const std::string &Cj  = Cs[j];
          const std::string &CjB = Cj+'B';
          if (b1.has_value(Cj)) {
            b1.setattr(CjB, Cj==keeper ? b1[Cj] : 0.0);
          }
        }

      } else if ( tb > 0 && ta == -2 ) {
        // just add the B state to the list of alchemical terms
        desres::msys::Destro &b1 = m1.block(alc_blk).append(m2[blk][tb]);
        b1.setattr( moiety, i->moiety );
	if ( b1.has_attr("ffio_ai") ) b1.setattr("ffio_ai", std::abs(i->ai));
	if ( b1.has_attr("ffio_aj") ) b1.setattr("ffio_aj", std::abs(i->aj));
	if ( b1.has_attr("ffio_ak") ) b1.setattr("ffio_ak", std::abs(i->ak));
	if ( b1.has_attr("ffio_al") ) b1.setattr("ffio_al", std::abs(i->al));
	if ( b1.has_attr("ffio_am") ) b1.setattr("ffio_am", std::abs(i->am));
	if ( b1.has_attr("ffio_an") ) b1.setattr("ffio_an", std::abs(i->an));
	if ( b1.has_attr("ffio_ao") ) b1.setattr("ffio_ao", std::abs(i->ao));
	if ( b1.has_attr("ffio_ap") ) b1.setattr("ffio_ap", std::abs(i->ap));
	for ( int j = 0; j < nCs; ++j ) {
          const std::string &Cj  = Cs[j];
          const std::string &CjB = Cj+'B';
          if (b1.has_value(Cj)) {
            b1.setattr(CjB, b1[Cj]);
          }
        }

      } else if ( tb > 0 && ta == -1 ) {
        // just add the B state to the list of regular terms.
        desres::msys::Destro &b1 = m1.block(blk).append( m2[blk][tb] );
	if ( b1.has_attr("ffio_ai") ) b1.setattr("ffio_ai", std::abs(i->ai));
	if ( b1.has_attr("ffio_aj") ) b1.setattr("ffio_aj", std::abs(i->aj));
	if ( b1.has_attr("ffio_ak") ) b1.setattr("ffio_ak", std::abs(i->ak));
	if ( b1.has_attr("ffio_al") ) b1.setattr("ffio_al", std::abs(i->al));
	if ( b1.has_attr("ffio_am") ) b1.setattr("ffio_am", std::abs(i->am));
	if ( b1.has_attr("ffio_an") ) b1.setattr("ffio_an", std::abs(i->an));
	if ( b1.has_attr("ffio_ao") ) b1.setattr("ffio_ao", std::abs(i->ao));
	if ( b1.has_attr("ffio_ap") ) b1.setattr("ffio_ap", std::abs(i->ap));
          
      } else if ( ta == 0 && tb > 0 ) {
        // A state has zeroed coefficients.  Need regular and alchemical term

        // build the A state
        desres::msys::Destro &b1 = m1.block(alc_blk).append(m2[blk][tb]);
        b1.setattr( moiety, i->moiety );
	if ( b1.has_attr("ffio_ai") ) b1.setattr("ffio_ai", std::abs(i->ai));
	if ( b1.has_attr("ffio_aj") ) b1.setattr("ffio_aj", std::abs(i->aj));
	if ( b1.has_attr("ffio_ak") ) b1.setattr("ffio_ak", std::abs(i->ak));
	if ( b1.has_attr("ffio_al") ) b1.setattr("ffio_al", std::abs(i->al));
	if ( b1.has_attr("ffio_am") ) b1.setattr("ffio_am", std::abs(i->am));
	if ( b1.has_attr("ffio_an") ) b1.setattr("ffio_an", std::abs(i->an));
	if ( b1.has_attr("ffio_ao") ) b1.setattr("ffio_ao", std::abs(i->ao));
	if ( b1.has_attr("ffio_ap") ) b1.setattr("ffio_ap", std::abs(i->ap));

	for ( int j = 0; j < nCs; ++j ) {
          const std::string &Cj  = Cs[j];
          const std::string &CjB = Cj+'B';
          if (b1.has_value(Cj)) {
            // build the B state
            b1.setattr(CjB, b1[Cj]);
            // zero the A state after building the A state
            if (Cj != keeper) b1.setattr(Cj, 0.0);
          }
        }

      } else if ( ta > 0 && tb > 0 ) {
        // move the A state to the alchemical block
        desres::msys::Destro &b1 = m1.block(alc_blk).append(m1[blk][ta]);
        b1.setattr( moiety, i->moiety );
        remove.insert(ta);
        // fetch data for the B state
        desres::msys::Destro const& b2 = m2[blk][tb];

        for (int j=0; j<nCs; ++j) {
          const std::string &Cj  = Cs[j];
          const std::string &CjB = Cj+'B';
          if (b2.has_value(Cj)) b1.setattr(CjB,b2[Cj]);
        }

      } else MSYS_FAIL("Unsupported mapping of " << blk << ", line ta(" << ta << ") and tb(" << tb << ")");
    }
    // remove the rows that were copied into the alchemical state
    for (std::set<int>::const_reverse_iterator i=remove.rbegin(); 
                                               i!=remove.rend(); ++i) {
      m1.block(blk).del(*i);
    }
  }

  /*! Fix atom indices for a constraint that was simply copied from stage 2
   * @param[out] cnew newly mapped constraint term
   * @param invmap map from original indices to alchemcally combined indices
   */
  void update_atom_indices( desres::msys::Destro &cnew,
                            std::map<int,int> &invmap ) {

    // search for elements of the form ffio_aX and change them to their
    // mapped value.
    typedef std::map<std::string,desres::msys::Destro::schema_t> SchemaMap;
    const SchemaMap &schemas = cnew.schemas();
    for (SchemaMap::const_iterator i=schemas.begin(); i!=schemas.end(); ++i) {
      std::string attr = i->first;
      if (!cnew.has_value(attr)) continue;
      if (!strncmp(attr.c_str(), "ffio_a", 6)) {
        cnew[attr] = invmap[cnew[attr]];
      }
    }
  }

  /*! Fix atom indices for constraint from stage 2 that was merged with stage 1
   * @param cold stage 1 constraint
   * @param cnew stage 2 constraint
   * @param invmap map from original indices to alchemcally combined indices
   */
  void append_atom_indices( desres::msys::Destro &cold, 
                            const desres::msys::Destro &cnew,
                            std::map<int,int> &invmap ) {

    // create a mapping from ffio_aj, ffio_ak, .... to ffio_c1, ffio_c2, ...
    static const char *from[] = {
      "ffio_aj", "ffio_ak", "ffio_al", "ffio_am", 
      "ffio_an", "ffio_ao", "ffio_ap", "ffio_aq"
    };
    static const char *to[] = {
      "ffio_c1", "ffio_c2", "ffio_c3", "ffio_c4", 
      "ffio_c5", "ffio_c6", "ffio_c7", "ffio_c8" 
    };
    static const unsigned a2c_size = sizeof(from)/sizeof(const char *);
    std::map<std::string,std::string> a2c;
    for (unsigned i=0; i<a2c_size; i++) a2c[from[i]] = to[i];

    // parse the ffio_funct records to get the number of light atoms.
    std::string oldfunct = CaseInsensitive::to_upper(cold("ffio_funct"));
    std::string newfunct = CaseInsensitive::to_upper(cnew("ffio_funct"));
    int oldnumH=0, newnumH=0;
    if (sscanf(oldfunct.c_str(), "AH%d", &oldnumH) != 1 ||
        sscanf(newfunct.c_str(), "AH%d", &newnumH) != 1)
      // nothing to do; we don't try to map these constraints
      return;
    if (oldnumH < 1 || newnumH < 1 ||
        oldnumH > (int)a2c_size || newnumH > (int)a2c_size) 
      MSYS_FAIL("Unsupported constraint type: " << oldfunct << ", " << newfunct);

    // cache the non-ai indices in a set
    std::set<int> indices;
    for (int i=0; i<oldnumH; i++) indices.insert(cold(from[i]));

    std::map<int,double> coeffs;
    // now go over the non-ai indices in the new record; if an atom in the
    // new record maps to an atom that's not in the old constraint, then
    // grab its coefficient.
    for (int i=0; i<newnumH; i++) {
      const std::string &attr = from[i];
      int ah = invmap[cnew[attr]];
      if (indices.find(ah) != indices.end()) {
        // this atom is already in the constraint
        continue;
      }
      // get the corresponding coefficient
      const std::string &coeff_attr = a2c[attr];
      coeffs[ah] = cnew(coeff_attr);
    }
    if (!coeffs.size()) return;
    if (indices.size() + coeffs.size() > a2c_size)
      MSYS_FAIL("Too many atoms in a constraint term");

    unsigned ind = indices.size();
    for (std::map<int,double>::const_iterator i=coeffs.begin();
                                              i!=coeffs.end(); ++i,++ind) {
      cold.setattr(from[ind], i->first);
      cold.setattr(to[ind], i->second);
    }
    // update ffio_funct
    char funct[32];
    int n = indices.size() + coeffs.size();
    sprintf(funct, "AH%d", n);
    cold.setattr("ffio_funct", funct);
  }

  void validate_ct( const desres::msys::Destro &m ) {
      using desres::msys::Destro;

      const Destro &atoms = m.block("m_atom");
      const Destro &ffio_ff = m.block("ffio_ff");
      const Destro &sites = ffio_ff.block("ffio_sites");

      int natoms=atoms.size();
      int nsites=sites.size();
      if (!nsites) MSYS_FAIL("Invalid ct: Empty ffio_sites block");
      int npseudos = 0;
      if (ffio_ff.has_block("ffio_pseudo")) {
          npseudos = ffio_ff.block("ffio_pseudo").size();
      }
      int nparticles=natoms+npseudos;
      int nblocks=nparticles/nsites;
      if (nblocks * nsites != nparticles) {
          MSYS_FAIL("Invalid ct: " << natoms << " atoms, " << npseudos << " pseudos, " << nsites << " sites");
      }
  }
  /*! combine vdw tables from stage 1 and 2.
   * @param m1 the stage 1 ct block
   * @param m2 the stage 2 ct block
   *
   * vdw types with names that 
   * don't appear in the stage 1 table are appended; vdw types with names
   * that do appear and have different values are given new names, and the
   * site that refer to those names are updated.
   */
  void vdw_combine( desres::msys::Destro &m1, desres::msys::Destro &m2 ) {

    desres::msys::Destro &fm1 = m1.block("ffio_ff");
    desres::msys::Destro &fm2 = m2.block("ffio_ff");

    if (!fm1.has_block( "ffio_vdwtypes" ) || 
        !fm2.has_block( "ffio_vdwtypes" ) ||
        !fm2.has_block( "ffio_sites" ) )
        return;

    desres::msys::Destro& vdw1 = fm1.block("ffio_vdwtypes");
    desres::msys::Destro const& vdw2 = fm2.block("ffio_vdwtypes");
    std::map<std::string,double> stage1_c1, stage1_c2;

    for (unsigned i=1; i<=vdw1.size(); i++) {
      const desres::msys::Destro &rec = vdw1[i];
      const std::string &name = rec("ffio_name");
      stage1_c1[name] = rec("ffio_c1").or_else(0.0);
      stage1_c2[name] = rec("ffio_c2").or_else(0.0);
    }
    for (unsigned i=1; i<=vdw2.size(); i++) {
      const desres::msys::Destro &rec = vdw2[i];
      const std::string &name = rec("ffio_name");
      if (stage1_c1.find( name ) == stage1_c1.end()) {
        // new vdw entry.
        vdw1.append( rec );
      } else {
        // compare c1, c2
        double c1 = rec("ffio_c1").or_else(0.0);
        double c2 = rec("ffio_c2").or_else(0.0);
        if (c1 != stage1_c1[name] || c2 != stage1_c2[name]) {
          // same name, different values.  Create a new entry and update sites
          desres::msys::Destro &new_entry = vdw1.append();
          std::string new_name = name + "_stage2";
          // be sure it's unique
          while ( stage1_c1.find( new_name ) != stage1_c1.end())
            new_name += "_";
          // stash values
          new_entry.setattr("ffio_name", new_name);
          new_entry.setattr("ffio_funct", rec("ffio_funct"));
          new_entry.setattr("ffio_c1", c1 );
          new_entry.setattr("ffio_c2", c2 );
          // update sites
          desres::msys::Destro &sites = fm2.block("ffio_sites");
          for (unsigned j=1; j<=sites.size(); j++) {
            desres::msys::Destro &site = sites[j];
            const std::string &type = site("ffio_vdwtype").or_else("");
            if (type == name) {
              site.setattr("ffio_vdwtype", new_name);
            }
          } // loop to update sites
        } // values differ 
      } // name conflict
    }
  }

  /*! Map m_bond records from stage2 into stage1
   *  @param a2inv mapping from ct2 atom index to 
   *  @param ct1 stage 1 ct
   *  @param ct2 stage 2 ct
   */
  void fixup_m_bond(const std::map<int,int>& a2inv, 
                    desres::msys::Destro &ct1,
                    const desres::msys::Destro &ct2) {

    // go through every m_bond record in ct2 and add it to ct1 after
    // mapping the atom indices.  
    if (!ct2.has_block("m_bond")) return; // nothing to do
    if (!ct1.has_block("m_bond")) ct1.new_array("m_bond");
    desres::msys::DestroArray &m_bond1 = ct1("m_bond");
    const desres::msys::DestroArray &m_bond2 = ct2("m_bond");

    // Keep track of the bonds we already have
    typedef std::set<std::pair<int,int> > BondSet;
    BondSet bondset;
    for (unsigned i=1; i<=m_bond1.size(); i++) {
      int from = m_bond1[i]("m_from");
      int to   = m_bond1[i]("m_to");
      bondset.insert(std::make_pair(from,to));
    }

    // add the others
    for (unsigned i=1; i<=m_bond2.size(); i++) {
      int from = m_bond2[i]("m_from");
      int to   = m_bond2[i]("m_to");
      std::map<int,int>::const_iterator from_iter = a2inv.find(from);
      std::map<int,int>::const_iterator to_iter   = a2inv.find(to);
      if (from_iter == a2inv.end() ||
          to_iter == a2inv.end()) 
        MSYS_FAIL("Missing entry in fepio_atommap for " << from << " " << to);

      BondSet::value_type p(from_iter->second, to_iter->second);
      if (bondset.find(p) != bondset.end()) continue;
      bondset.insert(p);
      desres::msys::Destro &elem = m_bond1.append();
      elem.setattr("m_from", p.first);
      elem.setattr("m_to",   p.second);
      elem.setattr("m_order",1);
    }
  }

  void combine_entry(desres::msys::Destro& m1, desres::msys::Destro const& m2,
                     std::string const& blk, std::map<int,int>& a2_inv_map) {

      // Make a hash of constraint terms in stage 1, indexed by the first atom.
      std::map<int,desres::msys::Destro *> stage1constraints;

      desres::msys::Destro & cons1 = 
        m1.has_block(blk) ? m1.block(blk) : m1.new_block(blk);

      for (unsigned i=0; i<cons1.size(); i++) {
        desres::msys::Destro &c = cons1[i+1];
        int ai = c("ffio_ai").or_else(0);
        if (ai<1) MSYS_FAIL("Invalid " << blk << " term found in stage 1");
        if (stage1constraints.find(ai) != stage1constraints.end())
          MSYS_FAIL("Found overlapping " << blk << " for atom " << ai);
        stage1constraints[ai] = &c;
      }
  
      desres::msys::Destro const& cons2 = m2.block(blk);
      for (unsigned i=0; i<cons2.size(); i++) {
        const desres::msys::Destro &c = cons2[i+1];
        int ai = c("ffio_ai").or_else(0);
        if (!ai) MSYS_FAIL("Invalid " << blk << " term in stage 2");
        int mapped_ai = a2_inv_map[ai];

        if (stage1constraints.find(mapped_ai) == stage1constraints.end()) {
          // new constraint term: copy it over and update the atom indices
          desres::msys::Destro &cnew = cons1.append(c);
          update_atom_indices( cnew, a2_inv_map );

        } else {
          // push stage2 atoms onto the existing constraint
          desres::msys::Destro &cold = *stage1constraints[mapped_ai];
          append_atom_indices( cold, c, a2_inv_map );
        }
      }
  }

  /*! append alchemical forms of atoms and groups onto mapped original forms
   * @param fm1 stage 1 ct
   * @param fm2 stage 2 ct
   * @param map fepio_fep map
   */
  void alchemical_combine( desres::msys::Destro &fm1, desres::msys::Destro const& fm2, 
                          const FepioMapping &map ) {
    desres::msys::Destro& atom_1 = fm1.block("m_atom");
    desres::msys::Destro const& atom_2 = fm2.block("m_atom");
    desres::msys::Destro& m1 = fm1.block("ffio_ff");
    desres::msys::Destro const& m2 = fm2.block("ffio_ff");

    FepioMapping::const_iterator atoms  = map.find("fepio_atommaps");
    FepioMapping::const_iterator bonds  = map.find("fepio_bondmaps");
    FepioMapping::const_iterator morses = map.find("fepio_morsebondmaps");
    FepioMapping::const_iterator angles = map.find("fepio_anglemaps");
    FepioMapping::const_iterator diheds = map.find("fepio_dihedmaps");
    FepioMapping::const_iterator pairs  = map.find("fepio_pairmaps");
    FepioMapping::const_iterator excls  = map.find("fepio_exclmaps");
    
    // Map perturbed ct atom numbers into combined ct numbers
    std::map<int,int> a2_inv_map;
    for (unsigned i=1; i<=atom_2.size(); i++) a2_inv_map[i]=i;

    const char *moiety = "ffio_moiety";
    if (atoms != map.end()) {
      for (FepList::const_iterator i=atoms->second.begin();
          i!=atoms->second.end(); ++i) {
        int const ai = i->ai;
        int const aj = i->aj;
        if (ai > 0 && aj > 0) {
          a2_inv_map[aj] = ai;
          desres::msys::Destro &t1 = m1["ffio_sites"][ai];
          t1.setattr( moiety, i->moiety );
          desres::msys::Destro const& t2 = m2["ffio_sites"][aj];
          t1.setattr("ffio_chargeB", double(t2["ffio_charge"]));
          t1.setattr("ffio_vdwtypeB", std::string(t2["ffio_vdwtype"]));
        } else if ( ai > 0 && aj < 0 ) {
          desres::msys::Destro &t1 = m1["ffio_sites"][ai];
          t1.setattr( moiety, i->moiety );
	  provide_dummy(m1);
	  t1.setattr("ffio_chargeB", 0.0);
	  t1.setattr("ffio_vdwtypeB", "DMY");
	} else if ( ai < 0 && aj > 0 ) {
	  a2_inv_map[aj] = std::abs(ai);
	  provide_dummy(m1);
          desres::msys::Destro const& b2 = m2["ffio_sites"][aj];
	  atom_1.append(atom_2[aj]);
	  desres::msys::Destro& ns1 = m1.block("ffio_sites").append(b2);
          ns1.setattr( moiety, i->moiety );
	  ns1.setattr("ffio_charge", 0.0);
	  ns1.setattr("ffio_chargeB", b2["ffio_charge"]);
	  ns1.setattr("ffio_vdwtype", "DMY");
	  ns1.setattr("ffio_vdwtypeB", b2["ffio_vdwtype"]);
	} else MSYS_FAIL("ai(" << ai << ") and aj(" << aj << ") < 0 in atommap");
      }
    }

    fixup_m_bond(a2_inv_map, fm1, fm2);

    if (bonds != map.end()) {
      const char *alc = "ffio_bonds_alchemical";
      m1.new_array(alc);
      combine_entry( m1, m2, bonds, "ffio_c1", "ffio_bonds", alc);
    }
    if (morses != map.end()) {
      const char *alc = "ffio_morsebonds_alchemical";
      m1.new_array(alc);
      combine_entry( m1, m2, morses, "ffio_c1", "ffio_morsebonds", alc);
    }
    if (angles != map.end()) {
      const char *alc = "ffio_angles_alchemical";
      m1.new_array(alc);
      combine_entry( m1, m2, angles, "", "ffio_angles", alc);
    }
    if (diheds != map.end()) {
      const char *alc = "ffio_dihedrals_alchemical";
      m1.new_array(alc);
      combine_entry( m1, m2, diheds, "", "ffio_dihedrals", alc);
    }
    if (pairs != map.end()) {
      const char *alc = "ffio_pairs_alchemical";
      m1.new_array(alc);
      combine_entry( m1, m2, pairs, "", "ffio_pairs", alc);
    }
    if (excls != map.end()) {
      unsigned row=1;
      for (FepList::const_iterator i=excls->second.begin();
          i!=excls->second.end(); ++i, ++row) {
        int ta = i->ti;
        int tb = i->tj;
        int ai = std::abs(i->ai);
        int aj = std::abs(i->aj);
        if (ai > aj) std::swap(ai, aj); // We require ai<aj in exclusion lists.
        else if (ai==aj) {
            MSYS_FAIL("Illegal self-exclusion at line " << row << " of fepio_exclmaps: ai=" 
                    << i->ai << " aj=" << i->aj);
        }
	if ( ta> 0 && tb == -1 ) {
	  // Nothing.
	} else if ( tb > 0 && ta == -1 ) {
	  desres::msys::Destro const& b2 = m2["ffio_exclusions"][tb];
	  desres::msys::Destro& b1 = m1.block("ffio_exclusions").append(b2);
	  b1.setattr("ffio_ai", ai);
	  b1.setattr("ffio_aj", aj);
	} else if ( ta == -1 && tb == -1 ) {
	  desres::msys::Destro& b1 = m1.block("ffio_exclusions").append();
	  b1.setattr("ffio_ai", ai);
	  b1.setattr("ffio_aj", aj);
	} else if ( ta > 0 && tb > 0 ) {
	  // nothing
	} else MSYS_FAIL("Unsupported exclusion line ta(" << ta << ") and tb(" << tb << ")");
      }
    }

    // constraints
    if (m2.has_block("ffio_constraints")) {
        combine_entry(m1, m2, "ffio_constraints", a2_inv_map);
    }
    // restraints handled the same way.  Note that we avoid creating 
    // an invalid restraint term with multiple atoms by virtual of the
    // fact that the funct of restraints does not start with "AH". 
    // It's a little fragile, but seems ok for now.
    // ported from viparr3,BRANCHES/3.5.5, DESRESCode#1380
    if (m2.has_block("ffio_restraints")) {
        combine_entry(m1, m2, "ffio_restraints", a2_inv_map);
    }

  }
}

  /*! perform alchemical combining if necessary
   * @param M_ original maeff file
   * @return copy of M_ if different from original
   */
namespace desres { namespace msys { namespace mae {
  std::string prep_alchemical_mae( const std::string& contents ) {
      size_t stage1=0, stage2=0;
      desres::msys::Maeff M(contents);
    
      for (unsigned i=0; i<M.size(); i++) {
        const desres::msys::Destro &ct = M[i+1];
        int stage = ct("fepio_stage").or_else(0);
        if (stage==1) stage1 = i+1;
        if (stage==2) stage2 = i+1;
      }
    
      if (stage1 && stage2) {
    
        // perform sanity checks on the alchemical blocks
        validate_ct(M[stage1]);
        validate_ct(M[stage2]);
    
        // combine vdw before touching sites since it might change vdwtype
        vdw_combine( M[stage1], M[stage2] );
    
        FepioMapping mapping;
        grovel_fepio_mapping( mapping, M[stage2]["fepio_fep"] );
        alchemical_combine( M[stage1], M[stage2], mapping );
        // FIXME - desres::msys::Maeff overrides del(std::string) but not del(size_t)
        desres::msys::Destro &m_ = M;
        m_.del( stage2 );
      }
      std::ostringstream out;
      out << M;
      return out.str();
    }
}}}

