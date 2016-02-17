/*
* Compiler problems with GCC on some systems...
* Doesn't like having cuserid called from C++
* So, this "C" stub was created to call it.
*
* Now it has also been rewritten to use current Posix standard
* for BSD compatibility and su -c compatibility
*/

#include <stdio.h>    
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

char* ccuserid (char* argument)
{
/*    return cuserid (argument); */

    struct passwd* pw = getpwuid ( geteuid() );
    return pw->pw_name;

}

