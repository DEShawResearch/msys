#ifndef WIN32
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

/* sockets */
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <errno.h>
#include <unistd.h>

static char* getpdb(const char* pdbcode) {
    struct hostent  *hp;
    struct sockaddr_in server;
    int hlen, rc, sock;
    char header[256];
    char filename[256];
    char buf[4096];

    static const char hostname[] = "files.rcsb.org";
    static const char user_agent[] = "msys/1.0";

    /* get host by name */
    if (!(hp=gethostbyname(hostname))) {
        fprintf(stderr, "Could not find PDB website: %s\n", hostname);
        return NULL;
    }

    char address[1030];
    sprintf(address, "%d.%d.%d.%d", 
            (unsigned char) hp->h_addr_list[0][0],
            (unsigned char) hp->h_addr_list[0][1],
            (unsigned char) hp->h_addr_list[0][2],
            (unsigned char) hp->h_addr_list[0][3]);
    memset(&server,0,sizeof(server));
    server.sin_family = PF_INET;
    server.sin_addr.s_addr = inet_addr(address);
    server.sin_port = htons(80);

    /* create socket */
    if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
        fprintf(stderr, "Error creating socket: %s\n", strerror(errno));
        return NULL;
    }

    /* connect */
    if (connect(sock, (struct sockaddr *)&server, sizeof(server)) < 0) {
        fprintf(stderr, "Error connecting to PDB website\n");
        return NULL;
    }

    /* send the header */
    sprintf(filename, "download/%s.pdb", pdbcode);
    sprintf(header, 
            "%s /%.256s HTTP/1.0\015\012User-Agent: %s\015\012%s\015\012", 
            "GET", filename, user_agent, "");
    hlen = strlen(header);

    rc = write(sock, header, hlen);
    if (rc!=hlen) {
        fprintf(stderr, "PDB write failed: %s\n", strerror(errno));
        close(sock);
        return NULL;
    }

    /* read everything returned.  */
    size_t nread=0, nmax=4096;
    char* data = malloc(nmax);
    while ((rc = read(sock, buf, sizeof(buf)))>0) {
        size_t newsize = nread + rc;
        if (newsize >= nmax) {
            nmax *= 2;
            data = realloc(data, nmax);
        }
        memcpy(data+nread, buf, rc);
        nread += rc;
    }
    data[nread]='\0';
    close(sock);
    return data;
}

static char* skip_header(char* data) {
    if (!data) return NULL;
    while (*data) {
        char* end = strchr(data, '\012');
        if (!end) return NULL;
        if (end-data==1) return end+1;
        data = end+1;
    }
    return NULL;
}

char* desres_msys_import_webpdb(const char* code) {
    char* data = getpdb(code);
    if (data) {
        char* pdb = skip_header(data);
        if (pdb) {
            pdb = strdup(pdb);
            free(data);
        }
        return pdb;
    }
    return NULL;
}
#endif
