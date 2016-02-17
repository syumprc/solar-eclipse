#include <stdlib.h>
#include <string.h>
#include <time.h>

void fdate_ (char* str, int stringlen)
{
    time_t now = time(0);
    char* timestr = ctime(&now);
    strncpy (str, timestr, (size_t) stringlen);
    int sindex = strlen(timestr);
    while (sindex < stringlen)
    {
	str[sindex++] = ' ';
    }
}
