#ifndef LIBMOLFILE_PLUGIN_H
#define LIBMOLFILE_PLUGIN_H
#include "vmdplugin.h"

#ifdef __cplusplus
extern "C" {
#endif
extern int molfile_dtrplugin_init(void);
extern int molfile_dtrplugin_register(void *, vmdplugin_register_cb);
extern int molfile_dtrplugin_fini(void);
extern int molfile_mol2plugin_init(void);
extern int molfile_mol2plugin_register(void *, vmdplugin_register_cb);
extern int molfile_mol2plugin_fini(void);
extern int molfile_parm7plugin_init(void);
extern int molfile_parm7plugin_register(void *, vmdplugin_register_cb);
extern int molfile_parm7plugin_fini(void);
extern int molfile_psfplugin_init(void);
extern int molfile_psfplugin_register(void *, vmdplugin_register_cb);
extern int molfile_psfplugin_fini(void);
extern int molfile_sdfplugin_init(void);
extern int molfile_sdfplugin_register(void *, vmdplugin_register_cb);
extern int molfile_sdfplugin_fini(void);
extern int molfile_dcdplugin_init(void);
extern int molfile_dcdplugin_register(void *, vmdplugin_register_cb);
extern int molfile_dcdplugin_fini(void);
extern int molfile_pdbplugin_init(void);
extern int molfile_pdbplugin_register(void *, vmdplugin_register_cb);
extern int molfile_pdbplugin_fini(void);
extern int molfile_webpdbplugin_init(void);
extern int molfile_webpdbplugin_register(void *, vmdplugin_register_cb);
extern int molfile_webpdbplugin_fini(void);
extern int molfile_xyzplugin_init(void);
extern int molfile_xyzplugin_register(void *, vmdplugin_register_cb);
extern int molfile_xyzplugin_fini(void);
extern int molfile_maeffplugin_init(void);
extern int molfile_maeffplugin_register(void *, vmdplugin_register_cb);
extern int molfile_maeffplugin_fini(void);
extern int molfile_dmsplugin_init(void);
extern int molfile_dmsplugin_register(void *, vmdplugin_register_cb);
extern int molfile_dmsplugin_fini(void);
#ifdef __cplusplus
}
#endif

#define MOLFILE_INIT_ALL \
molfile_dtrplugin_init(); \
molfile_mol2plugin_init(); \
molfile_parm7plugin_init(); \
molfile_psfplugin_init(); \
molfile_sdfplugin_init(); \
molfile_dcdplugin_init(); \
molfile_pdbplugin_init(); \
molfile_webpdbplugin_init(); \
molfile_xyzplugin_init(); \
molfile_maeffplugin_init(); \
molfile_dmsplugin_init(); \


#define MOLFILE_REGISTER_ALL(v,cb) \
molfile_dtrplugin_register(v,cb); \
molfile_mol2plugin_register(v,cb); \
molfile_parm7plugin_register(v,cb); \
molfile_psfplugin_register(v,cb); \
molfile_sdfplugin_register(v,cb); \
molfile_dcdplugin_register(v,cb); \
molfile_pdbplugin_register(v,cb); \
molfile_webpdbplugin_register(v,cb); \
molfile_xyzplugin_register(v,cb); \
molfile_maeffplugin_register(v,cb); \
molfile_dmsplugin_register(v,cb); \


#define MOLFILE_FINI_ALL \
molfile_dtrplugin_fini(); \
molfile_mol2plugin_fini(); \
molfile_parm7plugin_fini(); \
molfile_psfplugin_fini(); \
molfile_sdfplugin_fini(); \
molfile_dcdplugin_fini(); \
molfile_pdbplugin_fini(); \
molfile_webpdbplugin_fini(); \
molfile_xyzplugin_fini(); \
molfile_maeffplugin_fini(); \
molfile_dmsplugin_fini(); \


#endif
