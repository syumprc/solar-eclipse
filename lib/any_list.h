/* pipeshell.h
 * declarations for list functions
 * 
 *      Copyright (C) 1995 Southwest Foundation for Biomedical Research
 *                          All rights reserved.
 *                 Absolutely no warranty, express or implied.
 *
 * author:  Charles P. Peterson
 *   date:  June 8, 1995
 *
 */

#ifndef ANY_LIST_H
#define ANY_LIST_H

#include <stdlib.h>
#include <limits.h>
#include "safelib.h"

#define BOOLEAN int

#ifdef TRUE
#undef TRUE
#endif
#define TRUE 1

#ifdef FALSE
#undef FALSE
#endif
#define FALSE 0

struct any_list_st
{
    void *anyp;
    struct any_list_st *next;
};

typedef struct any_list_st any_list;
typedef any_list string_list;
typedef string_list String_List;

BOOLEAN any_list__add (any_list **list, void *anyp);
int any_list__count (any_list *list);
BOOLEAN any_list__remove (any_list **list, void *anyp);
BOOLEAN any_list__replace (any_list **list, void *anyp, void *newp);
void any_list__iterate_begin (any_list *list, any_list **temp);
BOOLEAN any_list__iterate (any_list **temp, void **anyp);
#define any_list__next any_list__iterate
#define any_list__begin any_list__iterate_begin
void any_list__delete_skeleton (any_list **list);
void *string_list__find (string_list *list, char *str);
/* 
 * Warning!  any_list__index is 1-based and inefficient!
 * Read notes in any_list.c before using!
 * Use any_list__begin and any_list__next for typical iteration.
 */
BOOLEAN any_list__index (any_list *list, int index, void **anyp);

#define List_Add any_list__add
#define List_Count any_list__count
#define List_RemoveFrom any_list__remove
#define List_Replace any_list__replace
#define List_Next any_list__iterate
#define List_Free any_list__delete_skeleton
#define List_Index any_list__index

#define String_List_Add any_list__add
#define String_List_Count any_list__count
#define String_List_RemoveFrom any_list__remove
#define String_List_Replace any_list__replace
#define String_List_Next any_list__iterate
#define String_List_Free any_list__delete_skeleton
#define String_List_Index any_list__index
#define String_List_Find string_list__find

#endif
