/*
 * expression.cc implements Expression and related classes
 *   and related classes
 * Written by Charles Peterson beginning on October 7, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#include <string.h>
#include <stdlib.h>
#include <math.h>        // pow
#include "safelib.h"
#include "any_list.h"
#include "token.h"
#include "expression.h"

Expression::Expression (const char *string)
{
    user_string = Strdup (string);

    Scan *s = new Scan (user_string);
    Token *nt = s->next_token ();
    Token *nt2 = s->next_token ();
    expr = nt->parse (nt2, NULL, s);
    delete s;
};

Expr *Binary_Expr::insert (Expr *new_expr, Expr **top)
{
    if (precedence () >= new_expr->precedence ())
    {
	 *top = new_expr;
	 return this;
    }
    else
    {
	*top = this;
	return second_expr->insert (new_expr, &second_expr);
    }
}

Named_Exp *Name_Or_Func_Expr (char *n, Scan *s)
{
    Token *nt = s->next_token ();
    Named_Exp *ne = 0;
    if (nt->type () == left_paren)
    {
	// This is a function call
	Token *start_args = new Begin_Token (s);
	Token *arg_tok = s->next_token ();
	Expr *arg = start_args->parse (arg_tok, NULL, s);
	Function_Expr *fe = new Function_Expr (n, arg);
	ne = fe;
    }
    else
    {
	// This is a variable; push token back
	s->push (nt);
	ne = new Name_Expr (n);
    }
    return ne;
}

Name_Node::Name_Node (const char *node_name, NAME_FUNCP node_fptr,
		      int node_key, Name_Node **name_list)
{
    name = Strdup (node_name);
    fptr = node_fptr;
    key = node_key;
    next = *name_list;
    *name_list = this;
}

void Name_Expr::bind (Context *c)
{
    name_function_id nfi;
    nfi = c->find (name);
    fptr = nfi.fptr;
    key = nfi.key;
}

name_function_id Context::find (char *name)
{
    Name_Node *search_list = name_list;
    while (search_list)
    {
	if (!StringCmp (search_list->name, name, false))  // Case insensitive
	{
	    name_function_id nfi;
	    nfi.fptr = search_list->fptr;
	    nfi.key = search_list->key;
	    return nfi;
	}
	search_list = search_list->next;
    }
    throw Undefined_Name (Strdup (name));
}


function_node *Function_Expr::function_list = NULL;
function_node *Function_Expr::it = NULL;

void function_setup (const char *name, int key, 
			   double (*fptr)(int*,double*,double*))
{
    function_node *fn = new function_node;
    fn->name = name;
    fn->key = key;
    fn->fptr = fptr;
    fn->next = Function_Expr::function_list;
    Function_Expr::function_list = fn;
}

Function_Expr::Function_Expr (char *n, Expr *argument) : Named_Exp (n)
{
    arg = argument;
    function_node *fn = function_list;
    while (fn)
    {
	if (!StringCmp (fn->name, n, case_ins))
	{
	    fptr = fn->fptr;
	    return;
	}
	fn = fn->next;
    }
    throw Undefined_Function();
}
