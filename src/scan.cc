/*
 * scan.cc implements Scan
 * This is the "scanner" used by token.cc
 * (A scanner breaks up a string into tokens used by the parser.)
 *
 * Written by Charles Peterson beginning on October 15, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "token.h"

#ifdef USE_SAFELIB
#include "safelib.h"
#else
#define Strdup strdup
#endif

Scan::Scan (char *string)
{
    scanbuf = string;
    ptr = scanbuf;
    token_list = NULL;
    begun = 0;
    paren_count = 0;
    all_tokens = 0;
}

Scan::~Scan ()
{
    while (all_tokens)
    {
	Token *next_token = all_tokens->next_all_tokens;
	delete all_tokens;
	all_tokens = next_token;
    }
}

Token *Scan::next_token (void)
{
    if (token_list)
    {
	Token *rt = token_list->token;
	Token_List *etl = token_list;
	token_list = token_list->next;
	delete etl;
	return rt;
    }

    if (!begun)
    {
	begun = 1;
	return new Begin_Token (this);
    }

    while (*ptr == ' ' || *ptr == '\t')
    { 
	ptr++; // Skip over whitespace (horiz)
    }

    char ch = *ptr;
    if (ch == '\0')
    {
	if (paren_count != 0) throw Syntax_Error();
	return new End_Token (this);
    }
    if (ch == '<')
    {
	if (ptr[1] == '=')
	{
	    ptr += 2;
	    return new Le (this);
	}
	if (ptr[1] == '<')
	{
	    ptr += 2;
	    return new Lt (this);
	}
	return scan_bracketed_name ();
    }

    if (ch == '>')
    {
	if (ptr[1] == '=')
	{
	    ptr += 2;
	    return new Ge (this);
	}
	if (ptr[1] == '>')
	{
	    ptr += 2;
	    return new Gt (this);
	}
    }

    if (ch == '=')
    {
	if (ptr[1] == '=')
	{
	    ptr += 2;
	    return new Eq (this);
	}
    }

    if (ch == '!')
    {
	if (ptr[1] == '=')
	{
	    ptr += 2;
	    return new Ne (this);
	}
    }
    if (ch == '.')
    {
	if (ptr[0] == '.' && tolower(ptr[1]) == 'e' &&
	    tolower(ptr[2]) == 'q' && ptr[3] == '.')
	{
	    ptr += 4;
	    return new Eq (this);
	}
	if (ptr[0] == '.' && tolower(ptr[1]) == 'n' &&
	    tolower(ptr[2]) == 'e' && ptr[3] == '.')
	{
	    ptr += 4;
	    return new Ne (this);
	}
	if (ptr[0] == '.' && tolower(ptr[1]) == 'g' &&
	    tolower(ptr[2]) == 'e' && ptr[3] == '.')
	{
	    ptr += 4;
	    return new Ge (this);
	}
	if (ptr[0] == '.' && tolower(ptr[1]) == 'l' &&
	    tolower(ptr[2]) == 'e' && ptr[3] == '.')
	{
	    ptr += 4;
	    return new Le (this);
	}
	if (ptr[0] == '.' && tolower(ptr[1]) == 'g' &&
	    tolower(ptr[2]) == 't' && ptr[3] == '.')
	{
	    ptr += 4;
	    return new Gt (this);
	}
	if (ptr[0] == '.' && tolower(ptr[1]) == 'l' &&
	    tolower(ptr[2]) == 't' && ptr[3] == '.')
	{
	    ptr += 4;
	    return new Lt (this);
	}
    }
    if (isdigit (ch) || ch == '.')
    {
	return scan_constant ();
    }
    if (isalpha (ch) || ch == '_')
    {
	return scan_name ();
    }
    if (ch == '*')
    {
	if (ptr[1] == '*')
	{
	    ptr += 2;
	    return new Double_Star (this);
	}
	ptr++;
	return new Star (this);
    }

// At this point, it must be a one-char token or an error

    ptr++;
    {
	switch (ch)
	{
	case '+':
	    return new Plus_Sign (this);
	case '-':
	    return new Minus_Sign (this);
	case '/':
	    return new Slash (this);
	case '^':
	    return new Hat (this);
	case '(':
	    paren_count++;
	    return new Left_Paren (this);
	case ')':
	    if (0 > --paren_count) throw Syntax_Error();
	    return new Right_Paren (this);
	}
    }
    throw Syntax_Error();
}

Token *Scan::scan_bracketed_name ()
{
    ptr++;
    char *eptr = strchr (ptr, '>');
    if (!eptr) throw Syntax_Error();
    char tch = *eptr;
    *eptr = '\0';
    Name *name = new Name (this);
    name->name = Strdup (ptr);
    *eptr = tch;
    ptr = eptr+1;
    return name;
}

Token *Scan::scan_name ()
{
    char *nptr = ptr;
    while (isalnum (*nptr) || *nptr == '_') nptr++;
    char tch = *nptr;
    *nptr = '\0';
    Name *name = new Name (this);
    name->name = Strdup (ptr);
    *nptr = tch;
    ptr = nptr;
    return name;
}

Token *Scan::scan_constant ()
{
// Allows leading + or -
// Allows exponent with E e D d (FORTRAN compatible) anywhere but beginning
// Allows . anywhere (even at beginning)

    char *cstring = NULL;

    { // Make a copy of the piece of string we need and advance ptr
	char *nptr = ptr;
	bool sign_allowed = false;  // leading - handled by parser
	bool exponent_allowed = true;
	bool decimal_allowed = true;
	while (isdigit (*nptr) || 
	       (sign_allowed && (*nptr == '-' || *nptr == '+')) ||
	       (exponent_allowed && (*nptr == 'e' || *nptr == 'E')) ||
	       (exponent_allowed && (*nptr == 'd' || *nptr == 'D')) ||
	       (decimal_allowed && (*nptr == '.')))
	{
	    if (*nptr == 'e' || *nptr == 'E' ||
		*nptr == 'd' || *nptr == 'D')
	    {
		decimal_allowed = false;
		exponent_allowed = false;
		sign_allowed = true;
	    }
	    else
	    {
		sign_allowed = false;
	    }
	    if (*nptr == '.')
	    {
		decimal_allowed = false;
	    }
	    nptr++;
	}
	char tch = *nptr;
	*nptr = '\0';
	cstring = Strdup (ptr);
	*nptr = tch;
	ptr = nptr;
    }
    
    char *cstr = cstring;
    while (*++cstr)
    {
	if (*cstr == 'd' || *cstr == 'D')
	{
	    *cstr = 'E';
	}
    }

    double value;
    char junk[1024];
    if (1 != sscanf (cstring, "%le%s", &value, &junk))
    {
	throw Syntax_Error();
    }

    free (cstring);
    Constant *con = new Constant (this, value);
    return con;
}

void Scan::push (Token *t)
{
    Token_List *tl = new Token_List;
    tl->next = token_list;
    tl->token = t;
    token_list = tl;
}
