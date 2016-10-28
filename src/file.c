/**
 * Copyright (c) 2012-2013, Mattias Frånberg
 * All rights reserved.
 *
 * This file is distributed under the Modified BSD License. See the COPYING file
 * for details.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <plinkio/file.h>

/**
 * Total size for a buffer containing the copy command
 * that will be issued to the shell.
 */
#define FILE_COPY_BUFFER_SIZE 4096

file_status_t
file_copy(const char *from_path, const char *to_path)
{
    char *copy_command = malloc( sizeof( char ) * ( strlen( from_path ) + strlen( to_path ) + 5 ) );
    int status = 0;

    sprintf( copy_command, "cp %s %s", from_path, to_path );
    status = system( copy_command );
    free( copy_command );

    if( status != -1 )
    {
        return FILE_OK;
    }
    else
    {
        return FILE_ERROR;
    }

}

file_status_t
file_remove(const char *path)
{
    char *rm_command = malloc( sizeof( char ) * ( strlen( path ) + 4 ) );
    int status = 0;
    
    sprintf( rm_command, "rm %s", path );
    status = system( rm_command );
    free( rm_command );

    if( status != -1 )
    {
        return FILE_OK;
    }
    else
    {
        return FILE_ERROR;
    }
}
