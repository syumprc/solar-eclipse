/*
 * token.h defines Token and related classes
 *
 * Written by Charles Peterson beginning on October 8, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#ifndef TOKEN_H
#define TOKEN_H

#include "expression.h"

enum ttype {begin_token, end_token, constant, name, left_paren, right_paren,
	    plus_sign, minus_sign, slash, star, double_star, hat, 
	    eq, ne, ge, le, gt, lt};

// Token(s) (defined delow)
class Token;

class Begin_Token;
class End_Token;

class Constant;
class Name;
class Left_Paren;
class Right_Paren;

class Binary_Op;
class Plus_Sign;
class Minus_Sign;
class Slash;
class Star;
class Double_Star;
class Hat;
class Eq;
class Ne;
class Ge;
class Le;
class Gt;
class Lt;

// Scanner
class Scan
{
    char begun;
    char *scanbuf;
    char *ptr;
    Token *scan_constant ();
    Token *scan_name ();
    Token *scan_bracketed_name ();
    int paren_count;

    class Token_List
    {
    public:
	Token *token;
	Token_List *next;
    };
    Token_List *token_list;

    friend class Token;
    Token *all_tokens;
public:
    Scan (char *string);
    ~Scan ();
    Token *next_token ();
    void push (Token *t);
};

class Token
{
    friend class Scan;
    Token *next_all_tokens;
public:
    Token (Scan *s) {next_all_tokens=s->all_tokens; s->all_tokens=this;}
    virtual ~Token () {}
    virtual ttype type () = 0;
    virtual Expr *parse (Token *t, Expr *e, Scan *s)
	{throw Syntax_Error();}  // Covers most cases!
};

class Constant : public Token
{
    double val;
public:
    Constant (Scan *s, double d = 0.0) : Token (s) {val = d;}
    ttype type () {return constant;}
    double value () {return val;}
// Uses default parse method.
// All sequences beginning with Constant are syntax errors.
};

class Name : public Token
{
public:
    Name (Scan *s) : Token (s) {}
    char *name;
    ttype type () {return ::name;}
// Uses default parse method.
// All sequences beginning with Constant are syntax errors.
};

class Left_Paren : public Token
{
public:
    Left_Paren (Scan *s) : Token (s) {}
    ttype type () {return left_paren;}
// Normally, '(' may only appear after Begin_Token or a binary operator
// EXCEPT if preceded by a Name, in which case it is a function call
// That determination is made by Name_Or_Func_Expr, which is a
// quasi virtual constructor for Name_Expr and Function_Expr.
};

class Right_Paren : public Token
{
public:
    Right_Paren (Scan *s) : Token (s) {}
    ttype type () {return right_paren;}
    Expr *parse (Token *t, Expr *e, Scan *s);
};

class Binary_Op : public Token
{
    Expr *parse_constant (Constant *t, Expr *e, Scan *s, bool negate=false);
    Expr *parse_name (Name *t, Expr *e, Scan *s, bool negate=false);
    Expr *parse_left_paren (Left_Paren *t, Expr *e, Scan *s, 
			    bool negate=false);
    Expr *parse_minus_sign (Minus_Sign *t, Expr *e, Scan *s);
    Expr *parse_plus_sign (Plus_Sign *t, Expr *e, Scan *s);
public:
    Binary_Op (Scan *s) : Token (s) {}
    virtual Binary_Expr *expr () = 0;  // Defined for each leaf class
    Expr *parse (Token *t, Expr *e, Scan *s) {
	switch (t->type ())
	{
	case constant:
	    return parse_constant ((Constant*) t, e, s);
	case name:
	    return parse_name ((Name*) t, e, s);
	case left_paren:
	    return parse_left_paren ((Left_Paren*) t, e, s);
	case plus_sign:
	    return parse_plus_sign ((Plus_Sign*) t, e, s);
	case minus_sign:
	    return parse_minus_sign ((Minus_Sign*) t, e, s);
	}
	throw Syntax_Error();
    }
};

class Plus_Sign : public Binary_Op
{
public:
    Plus_Sign (Scan *s) : Binary_Op (s) {}
    ttype type () {return plus_sign;}
    Binary_Expr *expr () {return new Addition_Expr;}
};

class Minus_Sign : public Binary_Op
{
public:
    Minus_Sign (Scan *s) : Binary_Op (s) {}
    ttype type () {return minus_sign;}
    Binary_Expr *expr () {return new Subtraction_Expr;}
};

class Star : public Binary_Op
{
public:
    Star (Scan *s) : Binary_Op (s) {}
    ttype type () {return star;}
    Binary_Expr *expr () {return new Multiplication_Expr;}
};

class Slash : public Binary_Op
{
public:
    Slash (Scan *s) : Binary_Op (s) {}
    ttype type () {return slash;}
    Binary_Expr *expr () {return new Division_Expr;}
};

class Double_Star : public Binary_Op
{
public: 
    Double_Star (Scan *s) : Binary_Op (s) {}
    ttype type () {return double_star;}
    Binary_Expr *expr () {return new Power_Expr;}
};

class Hat : public Binary_Op
{
public: 
    Hat (Scan *s) : Binary_Op (s) {}
    ttype type () {return hat;}
    Binary_Expr *expr () {return new Power_Expr;}
};

class Eq : public Binary_Op
{
public:
    Eq (Scan *s) : Binary_Op (s) {}
    ttype type () {return eq;}
    Binary_Expr *expr () {return new Eq_Expr;}
};

class Ne : public Binary_Op
{
public:
    Ne (Scan *s) : Binary_Op (s) {}
    ttype type () {return ne;}
    Binary_Expr *expr () {return new Ne_Expr;}
};

class Ge : public Binary_Op
{
public:
    Ge (Scan *s) : Binary_Op (s) {}
    ttype type () {return ge;}
    Binary_Expr *expr () {return new Ge_Expr;}
};

class Le : public Binary_Op
{
public:
    Le (Scan *s) : Binary_Op (s) {}
    ttype type () {return le;}
    Binary_Expr *expr () {return new Le_Expr;}
};

class Gt : public Binary_Op
{
public:
    Gt (Scan *s) : Binary_Op (s) {}
    ttype type () {return gt;}
    Binary_Expr *expr () {return new Gt_Expr;}
};

class Lt : public Binary_Op
{
public:
    Lt (Scan *s) : Binary_Op (s) {}
    ttype type () {return lt;}
    Binary_Expr *expr () {return new Lt_Expr;}
};

class Begin_Token : public Token
{
    Expr *parse_constant (Constant *t, Expr *e, Scan *s);
    Expr *parse_name (Name *t, Expr *e, Scan *s);
    Expr *parse_left_paren (Left_Paren *t, Expr *e, Scan *s);
    Expr *parse_plus_sign (Plus_Sign *t, Expr *e, Scan *s);
    Expr *parse_minus_sign (Minus_Sign *t, Expr *e, Scan *s);
public:
    Begin_Token (Scan *s) : Token (s) {}
    ttype type () {return begin_token;}
    virtual Expr *parse (Token *t, Expr *e, Scan *s) {
	switch (t->type ())
	{
	case constant:
	    return parse_constant ((Constant*) t, e, s);
	case name:
	    return parse_name ((Name*) t, e, s);
	case left_paren:
	    return parse_left_paren ((Left_Paren*) t, e, s);
	case plus_sign:
	    return parse_plus_sign ((Plus_Sign*) t, e, s);
	case minus_sign:
	    return parse_minus_sign ((Minus_Sign*) t, e, s);
	}
	throw Syntax_Error();
    }
};

class End_Token : public Token
{
public:
    End_Token (Scan *s) : Token (s) {}
    ttype type () {return end_token;}
    virtual Expr *parse (Token *t, Expr *e, Scan *s) {return e;}
};


#endif
