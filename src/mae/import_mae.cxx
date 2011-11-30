#include <boost/algorithm/string.hpp>
#include "mae.hxx"
#include "sitemap.hxx"
#include "vdwmap.hxx"
#include "ff.hxx"
#include "../mae.hxx"

#include "destro/prep_alchemical_mae.hxx"

#include <cstdio>
#include <fstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using desres::fastjson::Json;
using namespace desres::msys;
namespace bio = boost::iostreams;

namespace {

    static Id find_chain(SystemPtr h, const std::string& name) {
        IdList chains = h->chains();
        for (Id i=0, n=chains.size(); i<n; i++) {
            Id id = chains[i];
            if (h->chain(id).name==name) return id;
        }
        return BadId;
    }

    static Id find_residue(SystemPtr h, Id chn, const std::string& name,
                                                int resid ) {
        IdList residues = h->residuesForChain(chn);
        for (Id i=0, n=residues.size(); i<n; i++) {
            Id id = residues[i];
            if (h->residue(id).resid==resid && h->residue(id).name==name) {
                return id;
            }
        }
        return BadId;
    }

    const char * ANUMS     = "m_atomic_number";
    const char * RESIDS    = "m_residue_number";
    const char * RESNAMES  = "m_pdb_residue_name";
    const char * CHAINS    = "m_chain_name";
    const char * SEGIDS    = "m_pdb_segment_name";
    const char * NAMES     = "m_pdb_atom_name";
    const char * XCOL      = "m_x_coord";
    const char * YCOL      = "m_y_coord";
    const char * ZCOL      = "m_z_coord";
    const char * VXCOL     = "ffio_x_vel";
    const char * VYCOL     = "ffio_y_vel";
    const char * VZCOL     = "ffio_z_vel";
    const char * GRP_TEMP  = "ffio_grp_thermostat";
    const char * GRP_ENERGY = "ffio_grp_energy";
    const char * GRP_LIGAND = "ffio_grp_ligand";
    const char * GRP_BIAS   = "ffio_grp_cm_moi";
    const char * GRP_FROZEN = "ffio_grp_frozen";
    const char * FORMAL_CHG = "m_formal_charge";

    const std::string empty;

    bool contains_non_whitespace(const std::string& s) {
        for (unsigned i=0; i<s.size(); i++) if (!isspace(s[i])) return true;
        return false;
    }

    void import_cell( const Json& ct, SystemPtr h ) {
        h->global_cell.A[0]=ct.get("chorus_box_ax").as_float(0);
        h->global_cell.A[1]=ct.get("chorus_box_ay").as_float(0);
        h->global_cell.A[2]=ct.get("chorus_box_az").as_float(0);
        h->global_cell.B[0]=ct.get("chorus_box_bx").as_float(0);
        h->global_cell.B[1]=ct.get("chorus_box_by").as_float(0);
        h->global_cell.B[2]=ct.get("chorus_box_bz").as_float(0);
        h->global_cell.C[0]=ct.get("chorus_box_cx").as_float(0);
        h->global_cell.C[1]=ct.get("chorus_box_cy").as_float(0);
        h->global_cell.C[2]=ct.get("chorus_box_cz").as_float(0);
    }

    void import_particles( const Json& ct, SystemPtr h,
                           IdList& atoms,
                           int * natoms,
                           int * npseudos ) {
        *natoms=0;
        *npseudos=0;
        const Json& m_atom = ct.get("m_atom");
        if (!m_atom) return;

        Id chn = BadId;
        Id res = BadId;

        const Json& anums = m_atom.get(ANUMS);
        const Json& resids = m_atom.get(RESIDS);
        const Json& resnames = m_atom.get(RESNAMES);
        const Json& chains = m_atom.get(CHAINS);
        const Json& segids = m_atom.get(SEGIDS);
        const Json& names = m_atom.get(NAMES);
        const Json& x = m_atom.get(XCOL);
        const Json& y = m_atom.get(YCOL);
        const Json& z = m_atom.get(ZCOL);
        const Json& vx = m_atom.get(VXCOL);
        const Json& vy = m_atom.get(VYCOL);
        const Json& vz = m_atom.get(VZCOL);
        const Json& temp = m_atom.get(GRP_TEMP);
        const Json& nrg = m_atom.get(GRP_ENERGY);
        const Json& lig = m_atom.get(GRP_LIGAND);
        const Json& bias = m_atom.get(GRP_BIAS);
        const Json& frz = m_atom.get(GRP_FROZEN);
        const Json& formals = m_atom.get(FORMAL_CHG);
        Id gtmp=BadId, gene=BadId, glig=BadId, gbias=BadId, gfrz=BadId;
        Id seg=BadId;

        if (!!temp) gtmp=h->addAtomProp("grp_temperature", IntType);
        if (!!nrg)  gene=h->addAtomProp("grp_energy", IntType);
        if (!!lig)  glig=h->addAtomProp("grp_ligand", IntType);
        if (!!bias) gbias=h->addAtomProp("grp_bias", IntType);
        if (!!frz)  gfrz=h->addAtomProp("grp_frozen", IntType);
        if (!!segids) seg=h->addAtomProp("segid", StringType);

        int j,n = m_atom.get("__size__").as_int();
        for (j=0; j<n; j++) {
            int anum=anums.elem(j).as_int();
            std::string chainname=chains.elem(j).as_string("");
            int resid=resids.elem(j).as_int(0);
            std::string resname=resnames.elem(j).as_string("UNK");
            std::string name=names.elem(j).as_string("");

            boost::trim(chainname);
            boost::trim(resname);
            boost::trim(name);

            if (bad(chn) || h->chain(chn).name!=chainname) {
                chn = h->addChain();
                h->chain(chn).name = chainname;
                res = BadId;
            }

            if (bad(res) || h->residue(res).resid!=resid 
                         || h->residue(res).name!=resname) {
                res = h->addResidue(chn);
                h->residue(res).resid=resid;
                h->residue(res).name=resname;
            }
            Id id = h->addAtom(res);
            atom_t& atm = h->atom(id);
            atm.atomic_number = anum;
            atm.name = name;
            atm.x = x.elem(j).as_float(0);
            atm.y = y.elem(j).as_float(0);
            atm.z = z.elem(j).as_float(0);
            atm.vx = vx.elem(j).as_float(0);
            atm.vy = vy.elem(j).as_float(0);
            atm.vz = vz.elem(j).as_float(0);
            if (!!formals) atm.formal_charge = formals.elem(j).as_int(0);
            if (!!temp) h->atomPropValue(id,gtmp)=temp.elem(j).as_int(0);
            if (!!nrg)  h->atomPropValue(id,gene)=nrg.elem(j).as_int(0);
            if (!!lig)  h->atomPropValue(id,glig)=lig.elem(j).as_int(0);
            if (!!bias) h->atomPropValue(id,gbias)=bias.elem(j).as_int(0);
            if (!!frz)  h->atomPropValue(id,gfrz)=frz.elem(j).as_int(0);
            if (!!segids)  {
                std::string segid = segids.elem(j).as_string("");
                boost::trim(segid);
                h->atomPropValue(id,seg)=segid;
            }
            atoms.push_back(id);
            *natoms += 1;
        }

        const Json& m_bond = ct.get("m_bond");
        if (m_bond.valid()) {
            const Json& m_from = m_bond.get("m_from");
            const Json& m_to = m_bond.get("m_to");
            const Json& m_order = m_bond.get("m_order");
            const int natoms = atoms.size();
            int j,n = m_bond.get("__size__").as_int();
            for (j=0; j<n; j++) {
                int ai = m_from.elem(j).as_int();
                int aj = m_to.elem(j).as_int();
                if (ai<1 || aj<1 || ai>natoms || aj>natoms) {
                    std::ostringstream ss;
                    ss << "Bond " << j+1 << " between nonexistent atoms " 
                        << ai << ", " << aj;
                    throw std::runtime_error(ss.str());
                }
                if (ai==aj) {
                    std::ostringstream ss;
                    ss << "Bond " << j+1 << " to self " << ai;
                    throw std::runtime_error(ss.str());
                }
                if (ai>aj) {
                    int tmp=ai; ai=aj; aj=tmp;
                }
                Id a1 = atoms[ai-1];
                Id a2 = atoms[aj-1];
                Id bnd = h->addBond(a1,a2);
                h->bond(bnd).order = m_order.elem(j).as_int(1);
                if (h->bond(bnd).order==0) h->bond(bnd).order = 1;
            }
        }

        const Json& pseudo = ct.get("ffio_ff").get("ffio_pseudo");
        if (pseudo.valid()) {
            const Json& resids = pseudo.get("ffio_residue_number");
            const Json& resnames = pseudo.get("ffio_pdb_residue_name");
            const Json& chains = pseudo.get("ffio_chain_name");
            const Json& segids = pseudo.get("ffio_pdb_segment_name");
            const Json& names = pseudo.get("ffio_atom_name");
            const Json& x = pseudo.get("ffio_x_coord");
            const Json& y = pseudo.get("ffio_y_coord");
            const Json& z = pseudo.get("ffio_z_coord");
            const Json& vx = pseudo.get(VXCOL);
            const Json& vy = pseudo.get(VYCOL);
            const Json& vz = pseudo.get(VZCOL);

            chn = BadId;
            res = BadId;

            int j,n = pseudo.get("__size__").as_int();
            for (j=0; j<n; j++) {
                std::string segid = segids.elem(j).as_string("");
                std::string chainname=chains.elem(j).as_string("");
                int resid=resids.elem(j).as_int(0);
                std::string resname=resnames.elem(j).as_string("UNK");
                std::string name=names.elem(j).as_string("");

                boost::trim(segid);
                boost::trim(chainname);
                boost::trim(resname);
                boost::trim(name);

                if (bad(chn) || h->chain(chn).name!=chainname) {
                    chn = find_chain(h,chainname);
                    if (bad(chn)) {
                        chn = h->addChain();
                        h->chain(chn).name = chainname;
                    }
                    res = BadId;
                }

                if (bad(res) || h->residue(res).resid!=resid 
                             || h->residue(res).name!=resname) {
                    res = find_residue(h, chn, resname, resid);
                    if (bad(res)) {
                        res = h->addResidue(chn);
                        h->residue(res).resid= resid;
                        h->residue(res).name = resname;
                    }
                }

                Id id = h->addAtom(res);
                atom_t& atom = h->atom(id);
                atom.name = name;
                atom.x = x.elem(j).as_float(0);
                atom.y = y.elem(j).as_float(0);
                atom.z = z.elem(j).as_float(0);
                atom.vx = vx.elem(j).as_float(0);
                atom.vy = vy.elem(j).as_float(0);
                atom.vz = vz.elem(j).as_float(0);
                atoms.push_back(id);
                *npseudos += 1;
            }
        }
    }

    void write_ffinfo( const Json& ff, SystemPtr h ) {
        ParamTablePtr extra = h->auxTable("forcefield");
        if (!extra) {
            extra=ParamTable::create();
            extra->addProp( "id", IntType);
            extra->addProp( "path", StringType);
            extra->addProp( "info", StringType);
            h->addAuxTable("forcefield", extra);
        }
        std::string info;
        info += ff.get("viparr_command").as_string("");
        info += "\n";
        const Json& arr = ff.get("viparr_info").get("viparr_section");
        if (arr.valid()) for (int i=0; i<arr.size(); i++) {
            info += arr.elem(i).as_string("");
            info += "\n";
        }
        Id row = extra->addParam();
        extra->value(row,0) = row;
        extra->value(row,1) = ff.get("viparr_workdir").as_string("");
        extra->value(row,2) = info;
    }

    bool skippable(const std::string& s) {
        return s=="ffio_sites"
            || s=="ffio_atoms" 
            || s=="viparr_info"
            || s=="ffio_pseudo"
            ;
    }

    bool endswith(const char * word, const char * suffix ) {
        const char * a = strstr(word, suffix);
        if (!a) return false;
        size_t startlen = a-word;
        return startlen + strlen(suffix) == strlen(word);
    }
}
                           
namespace desres { namespace msys {

    SystemPtr ImportMAE( std::string const& path,
                         bool ignore_unrecognized ) {

        /* slurp in the file */
        std::ifstream file(path.c_str());
        if (!file) {
            std::stringstream ss;
            ss << "Failed opening MAE file at '" << path << "'";
            throw std::runtime_error(ss.str());
        }

        return ImportMAEFromStream(file, ignore_unrecognized);
    }

    SystemPtr ImportMAEFromStream( std::istream& file,
                                   bool ignore_unrecognized ) {

        bio::filtering_istream in;
        /* check for gzip magic number */
        if (file.get()==0x1f && file.get()==0x8b) {
            in.push(bio::gzip_decompressor());
        }
        file.seekg(0);
        in.push(file);

        std::stringstream buf;
        buf << in.rdbuf();

        /* create a json representation */
        Json M;
        mae::import_mae( buf, M );

        /* if alchemical, do the conversion on the original mae contents,
         * then recreate the json */
        int stage1=0, stage2=0;
        for (int i=0; i<M.size(); i++) {
            const Json& ct = M.elem(i);
            int stage = ct.get("fepio_stage").as_int(0);
            if (stage==1) stage1 = i+1;
            if (stage==2) stage2 = i+1;
        }
        if (stage1 && stage2) {
            std::string alc = prep_alchemical_mae(buf.str());
            std::istringstream in(alc);
            mae::import_mae( in, M );
        }

        SystemPtr h = System::create();

        for (int i=0; i<M.size(); i++) {
            const Json& ct = M.elem(i);
            const Json& type = ct.get("ffio_ct_type");
            if (type.valid() && !strcmp(type.as_string(), "full_system")) {
                continue;
            }
            
            IdList atoms;
            int natoms=0, npseudos=0;
            import_cell( ct, h );
            import_particles( ct, h, atoms, &natoms, &npseudos );

            const Json& ff = ct.get("ffio_ff");
            if (ff.valid()) {
                write_ffinfo(ff, h);

                const Json& sites = ff.get("ffio_atoms").valid() ?
                    ff.get("ffio_atoms") : ff.get("ffio_sites");

                mae::SiteMap sitemap( h, sites, atoms, natoms, npseudos );
                mae::VdwMap vdwmap( ff );

                int j,n = ff.size();
                for (j=0; j<n; j++) {
                    const Json& blk = ff.elem(j);
                    if (blk.kind()!=Json::Object) continue;
                    if (!blk.get("__size__").as_int()) continue;
                    std::string name = ff.key(j);
                    if (skippable(name)) continue;
                    const char * suffix = "_alchemical";
                    bool alchemical = endswith(name.c_str(), suffix);
                    if (alchemical) {
                        name = name.substr(0,name.size()-strlen(suffix));
                    }
                    const mae::Ffio * imp = mae::Ffio::get(name);
                    if (!imp) {
                        if (ignore_unrecognized) {
                            std::cerr << "skipping unrecognized block '" 
                                    << name << "'" << std::endl;
                        } else {
                            std::stringstream ss;
                            ss << "No handler for block '" << name << "'";
                            throw std::runtime_error(ss.str());
                        }
                    } else if (imp->wants_all()) {
                        imp->apply( h, ff,  sitemap, vdwmap, alchemical );
                    } else {
                        imp->apply( h, blk, sitemap, vdwmap, alchemical );
                    }
                }
            }
        }
        h->updateFragids();
        return h;
    }

}}
