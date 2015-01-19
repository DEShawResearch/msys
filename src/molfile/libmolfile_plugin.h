#ifndef LIBMOLFILE_PLUGIN_H
#define LIBMOLFILE_PLUGIN_H
#include "vmdplugin.h"

#ifdef __cplusplus
extern "C" {
#endif
int msys_plugin_register(void *, vmdplugin_register_cb);
int msys_dtrplugin_init(void);
int msys_dtrplugin_register(void *, vmdplugin_register_cb);
int msys_dcdplugin_init(void);
int msys_dcdplugin_register(void *, vmdplugin_register_cb);

#ifdef __cplusplus
}
#endif

#define MOLFILE_INIT_ALL \
msys_dtrplugin_init(); \
msys_dcdplugin_init(); \

#define MOLFILE_REGISTER_ALL(v,cb) \
msys_plugin_register(v,cb); \
msys_dtrplugin_register(v,cb); \
msys_dcdplugin_register(v,cb); \

#define MOLFILE_FINI_ALL \


#endif
