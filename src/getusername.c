#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <pwd.h>

int main()
{
    struct passwd* pw = getpwuid (geteuid());
    printf ("%s\n",pw->pw_name);
    return 0;
}
