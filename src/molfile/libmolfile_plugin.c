#include "libmolfile_plugin.h"

int msys_plugin_init(void);
int msys_plugin_register(void *, vmdplugin_register_cb);
int msys_dcdplugin_init(void);
int msys_dcdplugin_register(void *, vmdplugin_register_cb);
int msys_dtrplugin_init(void);
int msys_dtrplugin_register(void *, vmdplugin_register_cb);
int msys_dxplugin_init(void);
int msys_dxplugin_register(void *, vmdplugin_register_cb);
int msys_psfplugin_init(void);
int msys_psfplugin_register(void *, vmdplugin_register_cb);
int msys_xyzplugin_init(void);
int msys_xyzplugin_register(void *, vmdplugin_register_cb);

void molfile_init_all(void) {
    msys_plugin_init();
    msys_dcdplugin_init();
    msys_dtrplugin_init();
    msys_dxplugin_init();
    msys_psfplugin_init();
    msys_xyzplugin_init();
}

void molfile_register_all(void* v, vmdplugin_register_cb cb) {
    msys_plugin_register(v,cb);
    msys_dcdplugin_register(v,cb);
    msys_dtrplugin_register(v,cb);
    msys_dxplugin_register(v,cb);
    msys_psfplugin_register(v,cb);
    msys_xyzplugin_register(v,cb);
}

void molfile_fini_all(void) {
}

