
#ifndef CHARMM_FILE_H
#define CHARMM_FILE_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

int charmm_get_tokens(char **tok, int toklen,
			char *sbuf, int sbuflen,
			FILE *stream, int all_caps);

#ifdef __cplusplus
}
#endif

#endif

