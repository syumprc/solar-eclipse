/*
 * token.cc implements classes Token, Scannable, and Parsable
 * Written by Charles Peterson beginning on October 8, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#include <math.h>
#include "expression.h"
#include "token.h"

// Binary_Op methods

Expr *Binary_Op::parse_constant (Constant *c, Expr *e, Scan *s, bool negate)
{
    Binary_Expr *be = expr ();
    Expr *top = NULL;
    be->second_expr = new Constant_Expr (c->value ());
    be->first_expr = e->insert (be, &top);  // insert decides where top is

    if (negate)
    {
	Subtraction_Expr *se = new Subtraction_Expr;
	se->first_expr = new Constant_Expr (0.0);
	se->second_expr = be->second_expr;
	be->second_expr = se;
    }

    Token *nt = s->next_token ();
    Token *nt2 = s->next_token ();
    return nt->parse (nt2, top, s);
}

Expr *Binary_Op::parse_name (Name *t, Expr *e, Scan *s, bool negate)
{
    Binary_Expr *be = expr ();
    Expr *top = NULL;
    be->second_expr = Name_Or_Func_Expr (t->name, s);
    be->first_expr = e->insert (be, &top);  // insert decides where top is

    if (negate)
    {
	Subtraction_Expr *se = new Subtraction_Expr;
	se->first_expr = new Constant_Expr (0.0);
	se->second_expr = be->second_expr;
	be->second_expr = se;
    }

    Token *nt = s->next_token ();
    Token *nt2 = s->next_token ();
    return nt->parse (nt2, top, s);
}

Expr *Binary_Op::parse_left_paren (Left_Paren *t, Expr *e, Scan *s,
				   bool negate)
{
    Binary_Expr *be = expr ();

    Expr *top = NULL;
    be->first_expr = e->insert (be, &top);

    Token *bt = new Begin_Token (s);
    Token *nt = s->next_token ();
    Parenthesized_Expr *pe = new Parenthesized_Expr;
    be->second_expr = pe;
    pe->expr = bt->parse (nt, NULL, s);

    if (negate)
    {
	Subtraction_Expr *se = new Subtraction_Expr;
	se->first_expr = new Constant_Expr (0.0);
	se->second_expr = be->second_expr;
	be->second_expr = se;
    }

    Token *nt2x = s->next_token ();
    Token *nt3x = s->next_token ();
    
    return nt2x->parse (nt3x, top, s);
}

Expr *Binary_Op::parse_plus_sign (Plus_Sign *t, Expr *e, Scan *s)
{
// Ignore redundant "unary +" operator
// Simply get next token, and re-parse
    Token *nt = s->next_token ();
    return parse (nt, e, s);
}

Expr *Binary_Op::parse_minus_sign (Minus_Sign *t, Expr *e, Scan *s)
{
    Token *nt = s->next_token ();
    switch (nt->type ())
    {
    case constant:
	return parse_constant ((Constant *) nt, e, s, true);
    case name:
	return parse_name ((Name *) nt, e, s, true);
    case left_paren:
	return parse_left_paren ((Left_Paren *) nt, e, s, true);
    }
    throw Syntax_Error();
}

// Begin_Token methods

Expr *Begin_Token::parse_constant (Constant *c, Expr *e, Scan *s)
{
    Token *nt = s->next_token ();
    Token *nt2 = s->next_token ();

    return nt->parse (nt2, new Constant_Expr (c->value ()), s);
}

Expr *Begin_Token::parse_name (Name *n, Expr *e, Scan *s)
{
    Named_Exp *ne = Name_Or_Func_Expr (n->name, s);

    Token *nt = s->next_token ();
    Token *nt2 = s->next_token ();

    return nt->parse (nt2, ne, s);
}

Expr *Begin_Token::parse_left_paren (Left_Paren *lp, Expr *e, Scan *s)
{
    Token *nt = s->next_token ();
    Begin_Token *bt = new Begin_Token (s);
    Parenthesized_Expr *pe = new Parenthesized_Expr;
    pe->expr = bt->parse (nt, NULL, s);

    Token *nt2x = s->next_token ();
    Token *nt3x = s->next_token ();
    return nt2x->parse (nt3x, pe, s);
}

Expr *Begin_Token::parse_plus_sign (Plus_Sign *t, Expr *e, Scan *s)
{
// Just ignore it
    Token *nt = s->next_token ();
    return parse (nt, e, s);
}

Expr *Begin_Token::parse_minus_sign (Minus_Sign *t, Expr *e, Scan *s)
{
// Push 0 onto expr and evaluate as a "- <next token>" expression
    Expr *z = new Constant_Expr (0.0);
    Token *nt = s->next_token ();
    return t->parse (nt, z, s);
}

// Right_Paren methods

Expr *Right_Paren::parse (Token *t, Expr *e, Scan *s)
{
    s->push (t);
    return e;
}
