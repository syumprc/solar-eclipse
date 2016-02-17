#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <pwd.h>

int main()
{
    printf ("uid_t: %d\n", sizeof(uid_t));
    printf ("passwd: %d\n", sizeof (struct passwd));
    struct passwd* pw = getpwuid (geteuid());
    printf ("pw: %d\n", (int) pw);
    printf ("name: %s\n", pw->pw_name);

    return 0;
}
