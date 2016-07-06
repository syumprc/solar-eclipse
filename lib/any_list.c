/* pipeshell.h
 * list functions
 * 
 */

#include "any_list.h"

/* #define DEBUG__ANY_LIST */

BOOLEAN any_list__add (any_list **list, void *anyp)
{
    while (*list)
    {
	list = &(*list)->next;
    }
    *list = smalloc (sizeof (any_list));
    (*list)->anyp = anyp;
    (*list)->next = NULL;
    return 1;
}


int any_list__count (any_list *list)
{
    int count = 0;

    while (list)
    {
	list = list->next;
	count++;
    }
    return count;
}


BOOLEAN any_list__remove (any_list **list, void *anyp)
{
    while (*list)
    {
	if (anyp == (*list)->anyp)
	{
	    any_list *temp = *list;
	    *list = (*list)->next;
#ifdef DEBUG__ANY_LIST
    printf ("Removing node: %d   value: %d\n", 
	    (int) temp, (int) temp->anyp);
#endif
/*	    free (temp->anyp); not safe to do this yet */
	    free (temp);
	    return TRUE;
	}
	list = &(*list)->next;
    }
    return FALSE; /* No such value present */
}

BOOLEAN any_list__replace (any_list **list, void *anyp, void *newp)
{
    while (*list)
    {
	if (anyp == (*list)->anyp)
	{
	    (*list)->anyp = newp;
	    return TRUE;
	}
	list = &(*list)->next;
    }
    return FALSE; /* No such value present */
}

void any_list__iterate_begin (any_list *list, any_list **temp)
{
    *temp = list;
}

BOOLEAN any_list__iterate (any_list **temp, void **anyp)
{
    if (*temp)
    {
	*anyp = (*temp)->anyp;
	*temp = (*temp)->next;
	return TRUE;
    }
    else
    {
	*anyp = NULL;
	return FALSE;
    }
}

/*
 * NOTE: THIS USES A 1-BASED INDEX, NOT 0-BASED LIKE C
 *
 * Indexing is not a naturally efficient operation on lists
 * and we don't go out of our way to here to make it efficient.
 * Indexing lists should be used _very_ sparingly, such as when a
 * particular selection has been made from all listed items.
 *
 * iterate is the preferred operation, particularly when going through
 * the entire list.
 *
 * DO NOT USE index to select list members one after the other.  That is
 * what iterate is for.
 */
BOOLEAN any_list__index (any_list *list, int index, void **anyp)
{
    any_list *temp;
    int i;

    if (index < 1) return FALSE;

    any_list__iterate_begin (list, &temp);
    for (i = 1; i <= index; i++)
    {
	if (!any_list__iterate (&temp, anyp))
	{
	    return FALSE;
	}
    }
    return TRUE;
}

void any_list__delete_skeleton (any_list **list)
{
    while (*list)
    {
	any_list **this = list;
	list = &(*list)->next;
	free (*this);
    }
    *list = NULL;
}


void *string_list__find (string_list *list, char *str)
{
/*
 * The string list class adds the 'find' function, which finds a string
 * based on a matching name.  (Each other subclass of any_list could also
 * have distinct compare-by-value find.  A compare-by-identity (pointer
 * matching) find has not yet been needed, but it would be the obvious
 * the replace function--has not yet been needed, but would be the obvious
 * choice for the generic any_list.  But it doesn't seem that it would be
 * very useful, and it could easily be mistakenly used.
 *
 * Since it returns a pointer to the target string (or NULL), it could be
 * used in conjunction with either the replace or remove functions, which
 * take the pointer as input.
 */
    char *test_str = NULL;
    any_list *ilist = NULL;

    any_list__iterate_begin (list, &ilist);
    while (any_list__next (&ilist, (void **) &test_str))
    {
	if (!strcmp (test_str, str))
	{
	    return (char *) test_str;
	}
    }
    return NULL;
}
