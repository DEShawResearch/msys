#ifndef LIBMOLFILE_PLUGIN_H
#define LIBMOLFILE_PLUGIN_H
#include "vmdplugin.h"

#ifdef __cplusplus
extern "C" {
#endif
extern int msys_plugin_register(void *, vmdplugin_register_cb);
#ifdef __cplusplus
}
#endif

#define MOLFILE_INIT_ALL \


#define MOLFILE_REGISTER_ALL(v,cb) \
msys_plugin_register(v,cb); \


#define MOLFILE_FINI_ALL \


#endif
