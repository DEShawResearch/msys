#ifndef LIBMOLFILE_PLUGIN_H
#define LIBMOLFILE_PLUGIN_H

#include "vmdplugin.h"

#ifdef __cplusplus
extern "C" {
#endif
void molfile_init_all(void);
void molfile_register_all(void*, vmdplugin_register_cb);
void molfile_fini_all(void);
#ifdef __cplusplus
}
#endif

/* Backwards-compatibility macros. */
#define MOLFILE_INIT_ALL            molfile_init_all();
#define MOLFILE_REGISTER_ALL(v,cb)  molfile_register_all(v,cb);
#define MOLFILE_FINI_ALL            molfile_fini_all();

#endif
