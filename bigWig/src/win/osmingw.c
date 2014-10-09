#include "common.h"
#include <io.h>
#include <direct.h>
#include <dirent.h>
#include <sys/utime.h>
#undef BYTE
#undef WORD
#undef UWORD
#undef boolean
#include <winsock2.h>
#include <ws2tcpip.h>
#include "portable.h"

void sleep1000(int milli)
/* Sleep for given number of 1000ths of second */
{
if (milli > 0)
    {
    struct timeval tv;
    tv.tv_sec = milli/1000;
    tv.tv_usec = (milli%1000)*1000;
    select(0, NULL, NULL, NULL, &tv);
    }
}

struct fileInfo *newFileInfo(char *name, off_t size, bool isDir, int statErrno, 
	time_t lastAccess)
/* Return a new fileInfo. */
{
int len = strlen(name);
struct fileInfo *fi = needMem(sizeof(*fi) + len);
fi->size = size;
fi->isDir = isDir;
fi->statErrno = statErrno;
fi->lastAccess = lastAccess;
strcpy(fi->name, name);
return fi;
}

int cmpFileInfo(const void *va, const void *vb)
/* Compare two fileInfo. */
{
const struct fileInfo *a = *((struct fileInfo **)va);
const struct fileInfo *b = *((struct fileInfo **)vb);
return strcmp(a->name, b->name);
}

boolean makeDir(char *dirName)
/* Make dir.  Returns TRUE on success.  Returns FALSE
 * if failed because directory exists.  Prints error
 * message and aborts on other error. */
{
int err;
if ((err = _mkdir(dirName)) < 0)
    {
    if (errno != EEXIST)
	{
	perror("");
	errAbort("Couldn't make directory %s", dirName);
	}
    return FALSE;
    }
return TRUE;
}


struct fileInfo *listDirXExt(char *dir, char *pattern, boolean fullPath, boolean ignoreStatFailures)
/* Return list of files matching wildcard pattern with
 * extra info. If full path is true then the path will be
 * included in the name of each file. */
{
struct fileInfo *list = NULL, *el;
struct dirent *de;
DIR *d;
int dirNameSize = strlen(dir);
int fileNameOffset = dirNameSize+1;
char pathName[512];

if ((d = opendir(dir)) == NULL)
    return NULL;
memcpy(pathName, dir, dirNameSize);
pathName[dirNameSize] = '/';

while ((de = readdir(d)) != NULL)
    {
    char *fileName = de->d_name;
    if (differentString(fileName, ".") && differentString(fileName, ".."))
	{
	if (pattern == NULL || wildMatch(pattern, fileName))
	    {
	    struct stat st;
	    bool isDir = FALSE;
	    int statErrno = 0;
	    strcpy(pathName+fileNameOffset, fileName);
	    if (stat(pathName, &st) < 0)
		{
		if (ignoreStatFailures)
		    statErrno = errno;
		else
    		    errAbort("stat failed in listDirX");
		}
	    if (S_ISDIR(st.st_mode))
		isDir = TRUE;
	    if (fullPath)
		fileName = pathName;
	    el = newFileInfo(fileName, st.st_size, isDir, statErrno, st.st_atime);
	    slAddHead(&list, el);
	    }
	}
    }
closedir(d);
slSort(&list, cmpFileInfo);
return list;
}

struct fileInfo *listDirX(char *dir, char *pattern, boolean fullPath)
/* Return list of files matching wildcard pattern with
 * extra info. If full path is true then the path will be
 * included in the name of each file. */
{
return listDirXExt(dir, pattern, fullPath, FALSE);
}

time_t fileModTime(char *pathName)
/* Return file last modification time.  The units of
 * these may vary from OS to OS, but you can depend on
 * later files having a larger time. */
{
struct stat st;
if (stat(pathName, &st) < 0)
    errAbort("stat failed in fileModTime: %s", pathName);
return st.st_mtime;
}

char *getHost()
/* Return host name. */
{
static char *hostName = NULL;
static char buf[128];
if (hostName == NULL)
    {
    gethostname(buf, 127);
    strncpy(buf, hostName, sizeof(buf));
    chopSuffix(buf);
    hostName = buf;
    }
return hostName;
}

char *semiUniqName(char *base)
/* Figure out a name likely to be unique.
 * Name will have no periods.  Returns a static
 * buffer, so best to clone result unless using
 * immediately. */
{
int pid = getpid();
int num = time(NULL)&0xFFFFF;
char host[512];
strcpy(host, getHost());
char *s = strchr(host, '.');
if (s != NULL)
     *s = 0;
subChar(host, '-', '_');
subChar(host, ':', '_');
static char name[PATH_LEN];
safef(name, sizeof(name), "%s_%s_%x_%x",
	base, host, pid, num);
return name;
}

char *rTempName(char *dir, char *base, char *suffix)
/* Make a temp name that's almost certainly unique. */
{
char *x;
static char fileName[PATH_LEN];
int i;
char *lastSlash = (lastChar(dir) == '/' ? "" : "/");
for (i=0;;++i)
    {
    x = semiUniqName(base);
    safef(fileName, sizeof(fileName), "%s%s%s%d%s",
    	dir, lastSlash, x, i, suffix);
    if (!fileExists(fileName))
        break;
    }
return fileName;
}


void childExecFailedExit(char *msg)
/* Child exec failed, so quit without atexit cleanup */
{
fprintf(stderr, "child exec failed: %s\n", msg);
fflush(stderr);
_exit(1);  // Let the parent know that the child failed by returning 1.

/* Explanation:
_exit() is not the normal exit().  
_exit() avoids the usual atexit() cleanup.
The MySQL library that we link to uses atexit() cleanup to close any open MySql connections.
However, because the child's mysql connections are shared by the parent,
this causes the parent MySQL connections to become invalid,
and causes the puzzling "MySQL has gone away" error in the parent
when it tries to use its now invalid MySQL connections.
*/

}

void vaDumpStack(char *format, va_list args)
/* debugging function to run the pstack program on the current process. In
 * prints a message, following by a new line, and then the stack track.  Just
 * prints errors to stderr rather than aborts. For debugging purposes
 * only.  */
{
	fprintf(stderr, "dump stack not implemented on MinGW\n");
}

void dumpStack(char *format, ...)
/* debugging function to run the pstack program on the current process. In
 * prints a message, following by a new line, and then the stack track.  Just
 * prints errors to stderr rather than aborts. For debugging purposes
 * only.  */
{
va_list args;
va_start(args, format);
vaDumpStack(format, args);
va_end(args);
}

boolean maybeTouchFile(char *fileName)
/* If file exists, set its access and mod times to now.  If it doesn't exist, create it.
 * Return FALSE if we have a problem doing so (e.g. when qateam is gdb'ing and code tries 
 * to touch some file owned by www). */
{
if (fileExists(fileName))
    {
    struct utimbuf ut;
    ut.actime = ut.modtime = clock1();
    int ret = utime(fileName, &ut);
    if (ret != 0)
	{
	warn("utime(%s) failed (ownership?)", fileName);
	return FALSE;
	}
    }
else
    {
    FILE *f = fopen(fileName, "w");
    if (f == NULL)
	return FALSE;
    else
	carefulClose(&f);
    }
return TRUE;
}


/* Implement setenv, unsetenv since they don't exist on windows */
int setenv(const char *name, const char *value, int overwrite) {
	int len;
	char *str;

	if (overwrite == 0) {
		char * var;

		var = getenv(name);

		if (var != NULL)
			return 0;
	}

	len = strlen(name)+1+strlen(value)+1; 
    str = malloc(len*sizeof(char)); 
    sprintf(str, "%s=%s", name, value); 
    if(putenv(str) != EXIT_SUCCESS) {
      free(str);
      return -1;
    }

  	return 0;
}

int unsetenv(const char *name) {
	int len = strlen(name)+1+1;
	char * str = malloc(len*sizeof(char));

	/* On MINGW, a call `putenv("FOO=")' removes `FOO' from the 
     environment, rather than inserting it with an empty value. 
	*/
	sprintf(str, "%s=", name);

	if (putenv(str) != EXIT_SUCCESS) {
		free(str);
		return -1;
	}
	free(str);
	return 0;
}


/* pipe for net.c function */
int pipe(int pipefd[2]) {
	return _pipe( pipefd, 65536, O_BINARY ); /* Same capacity as Linux since 2.6.11 */ 
}
