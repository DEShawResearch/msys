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
    const char * INSERTION =  "m_insertion_code";
    const char * GROW_NAME =  "m_grow_name";

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
        const Json& inserts = m_atom.get(INSERTION);
        const Json& grow = m_atom.get(GROW_NAME);

        Id gtmp=BadId, gene=BadId, glig=BadId, gbias=BadId, gfrz=BadId;
        Id growcol=BadId;

        if (!!temp) gtmp=h->addAtomProp("grp_temperature", IntType);
        if (!!nrg)  gene=h->addAtomProp("grp_energy", IntType);
        if (!!lig)  glig=h->addAtomProp("grp_ligand", IntType);
        if (!!bias) gbias=h->addAtomProp("grp_bias", IntType);
        if (!!frz)  gfrz=h->addAtomProp("grp_frozen", IntType);

        SystemImporter imp(h);

        int j,n = m_atom.get("__size__").as_int();
        for (j=0; j<n; j++) {
            int anum=anums.elem(j).as_int();
            std::string chainname=chains.elem(j).as_string("");
            int resid=resids.elem(j).as_int(0);
            std::string resname=resnames.elem(j).as_string("UNK");
            std::string name=names.elem(j).as_string("");
            std::string segid=segids.elem(j).as_string("");
            const char* insert=inserts.elem(j).as_string("");

            Id id = imp.addAtom(chainname, segid, resid, resname, name, insert);

            atom_t& atm = h->atom(id);
            atm.atomic_number = anum;
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
            if (!!grow) {
                std::string g(grow.elem(j).as_string(""));
                boost::trim(g);
                if (!g.empty()) {
                    if (bad(growcol)) {
                        growcol=h->addAtomProp("m_grow_name", StringType);
                    }
                    h->atomPropValue(id,growcol)=g;
                }
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
                h->bond(bnd).resonant_order = h->bond(bnd).order;
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

            int j,n = pseudo.get("__size__").as_int();
            for (j=0; j<n; j++) {
                std::string segid = segids.elem(j).as_string("");
                std::string chainname=chains.elem(j).as_string("");
                int resid=resids.elem(j).as_int(0);
                std::string resname=resnames.elem(j).as_string("UNK");
                std::string name=names.elem(j).as_string("");

                Id id = imp.addAtom(chainname, segid, resid, resname, name);
                atom_t& atom = h->atom(id);
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
        /* hash by path so we don't add duplicates */
        std::map<std::string,Id> pathmap;
        for (Id i=0; i<extra->paramCount(); i++) {
            pathmap[extra->value(i,1).asString()] = i;
        }
        /* parse the viparr command to get the forcefield paths */
        std::string cmd = ff.get("viparr_command").as_string("");
        std::string workdir = ff.get("viparr_workdir").as_string(".");
        workdir += "/";
        std::vector<std::string> tokens;
        boost::split(tokens, cmd, boost::is_any_of(" "));
        if (tokens.size()>2) {
            for (unsigned i=1; i<tokens.size()-1; i++) {
                if (tokens[i]=="-f") {
                    std::string const& path = tokens.at(i+1);
                    if (!pathmap.count(path)) {
                        Id row = extra->addParam();
                        extra->value(row,0) = row;
                        extra->value(row,1) = path;
                        pathmap[path] = row;
                    }
                } else if (tokens[i]=="-d" || tokens[i]=="-m") {
                    std::string path = tokens.at(i+1);
                    if (path.substr(0,1)!="/") {
                        path = workdir + path;
                    }
                    if (!pathmap.count(path)) {
                        Id row = extra->addParam();
                        extra->value(row,0) = row;
                        extra->value(row,1) = path;
                        pathmap[path] = row;
                    }
                }
            }
        }
        /* add a provenance entry */
        if (h->provenance().empty()) {
            Provenance p;
            p.timestamp = ff.get("date").as_string("");
            p.user = ff.get("user").as_string("");
            p.workdir = ff.get("viparr_workdir").as_string("");
            p.cmdline = ff.get("viparr_command").as_string("");
            h->addProvenance(p);
        }
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

    bool is_full_system(Json const& ct) {
        const Json& type = ct.get("ffio_ct_type");
        return type.valid() && !strcmp(type.as_string(), "full_system");
    }

    void append_system(SystemPtr h, Json const& ct,
                       const bool ignore_unrecognized,
                       const bool structure_only) {

        h->name = ct.get("m_title").as_string("");

        IdList atoms;
        int natoms=0, npseudos=0;
        import_cell( ct, h );
        import_particles( ct, h, atoms, &natoms, &npseudos );
        if (structure_only) return;

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
                    imp->apply( h, ff,  sitemap, vdwmap );
                } else {
                    imp->apply( h, blk, sitemap, vdwmap );
                }
            }
        }
    }

    class iterator : public LoadIterator {
        std::stringstream buf;
        Json M;
        int ctnumber;
        const bool ignore_unrecognized;
        const bool structure_only;

    public:
        iterator(std::istream& file,
                 bool _ignore_unrecognized, bool _structure_only)
        : ctnumber(), 
          ignore_unrecognized(_ignore_unrecognized), 
          structure_only(_structure_only) {

            bio::filtering_istream in;
            /* check for gzip magic number */
            if (file.get()==0x1f && file.get()==0x8b) {
                in.push(bio::gzip_decompressor());
            }
            file.seekg(0);
            in.push(file);

            buf << in.rdbuf();

            /* create a json representation */
            mae::import_mae( buf, M );
            ctnumber=0;
        }

        SystemPtr next();
        SystemPtr read_all();
    };

    SystemPtr iterator::next() {
        for (; ctnumber<M.size(); ++ctnumber) {
            if (is_full_system(M.elem(ctnumber))) continue;
            SystemPtr h = System::create();
            append_system(h, M.elem(ctnumber++), 
                    ignore_unrecognized, structure_only);
            h->analyze();
            return h;
        }
        return SystemPtr();
    }

    SystemPtr iterator::read_all() {

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

        /* read all ct blocks into the same system */
        SystemPtr h = System::create();
        for (int i=0; i<M.size(); i++) {
            const Json& ct = M.elem(i);
            if (is_full_system(ct)) continue;
            append_system(h, ct, ignore_unrecognized, structure_only);
        }
        h->analyze();
        return h;
    }
}
                           
namespace desres { namespace msys {

    SystemPtr ImportMAE( std::string const& path,
                         bool ignore_unrecognized,
                         bool structure_only ) {

        std::ifstream file(path.c_str());
        if (!file) {
            MSYS_FAIL("Failed opening MAE file at '" << path << "'");
        }
        iterator it(file, ignore_unrecognized, structure_only);
        SystemPtr sys = it.read_all();
        sys->name = path;
        return sys;
    }

    SystemPtr ImportMAEFromBytes( const char* bytes, int64_t len,
                         bool ignore_unrecognized, bool structure_only ) {

        std::istringstream file;
        file.rdbuf()->pubsetbuf(const_cast<char *>(bytes), len);
        return iterator (file, ignore_unrecognized, structure_only).read_all();
    }

    SystemPtr ImportMAEFromStream( std::istream& file,
                                   bool ignore_unrecognized,
                                   bool structure_only) {

        return iterator (file, ignore_unrecognized, structure_only).read_all();
    }

    LoadIteratorPtr MaeIterator(std::string const& path) {
        const bool ignore_unrecognized = false;
        const bool structure_only = false;
        std::ifstream file(path.c_str());
        if (!file) {
            MSYS_FAIL("Failed opening MAE file at '" << path << "'");
        }
        return LoadIteratorPtr(
                new iterator(file, ignore_unrecognized, structure_only));
    }

}}
