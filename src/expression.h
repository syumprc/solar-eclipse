/*
 * expression.h defines class Expression:
 *    a interface class for expression parsing and evaluation
 * Related classes are declared and/or defined in:
 *    token.h expression.cc token.cc
 *
 * Written by Charles Peterson beginning on October 8, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#ifndef EXPRESSION_H
#define EXPRESSION_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

// bool type not defined in some compilers
// but it's a reserved word (not a definition) in ANSI C++ compilers
// Makefile must define NEEDS_BOOL for pre-ANSI compilers
#ifdef NEEDS_BOOL
enum bool {false, true};
#undef NEEDS_BOOL
#endif
#define Istrue(iexp)(((iexp)==0)?false:true)


// Exception classes
class Syntax_Error {};
class Undefined_Function {}; 
class Unresolved_Name {};
class Undefined_Name 
{
public:
    char *name;
    Undefined_Name (char *n) {name = n;}
};

// Precedences from lowest to highest
enum {Eq_Precedence, Add_Precedence, Mult_Precedence, Power_Precedence, 
      Default_Precedence};

typedef double (*NAME_FUNCP) (int key);

struct name_function_id
{
    NAME_FUNCP fptr;
    int key;
};

class Name_Node  // Element in a (name) Context
{
    char *name;
    NAME_FUNCP fptr;
    int key;
    Name_Node *next;

    Name_Node (const char *name, NAME_FUNCP fptr, int key,
	       Name_Node **name_list);
    ~Name_Node () {free (name);}
    friend class Context;
};

class Context
{
    Name_Node *name_list;
public:
    Context () {name_list = 0;}
    ~Context () {while (name_list)
    {Name_Node *n = name_list->next; delete name_list; name_list = n;}}
    void add (const char *name, NAME_FUNCP fptr, int key=0)
	{new Name_Node (name, fptr, key, &name_list);}
    name_function_id find (char *name);
};

class Expr
{
public:
    virtual ~Expr () {};
    virtual int precedence () {return Default_Precedence;}
    friend class Expression;
    static Expr *build (char *string);
    virtual double eval () = 0;
    virtual int is_empty () {return 0;}
    virtual Expr *insert (Expr *new_expr, Expr **top)
	{*top = new_expr; return this;}  // Default precedence handling
					 // (left to right)
    virtual void bind (Context *c) = 0;
    virtual char *debug (char *buf) = 0;
    virtual bool boolean () {return false;}
};

class Expression  // User wrapper around Expr
{
    Expr *expr;
    char *user_string;
public:
    Expression (const char *string);
    ~Expression () {delete user_string; delete expr;}
    char *string () {return user_string;}
    double eval () {return expr->eval();}
    char *debug (char *buf) {return expr->debug (buf);}
    bool boolean() {return expr->boolean();}
    void bind (Context *c) {expr->bind (c);}
};

class Empty_Expr : public Expr
{
public:
    virtual int is_empty () {return 1;}
    double eval () {return 0.0;} // This allows a leading unary minus
    void bind (Context *c) {}					
    char *debug (char *buf) {return strcpy (buf,"()");}
};

class Constant_Expr : public Expr
{
    double value;
public:
    Constant_Expr (double d) {value=d;}
    double eval () {return value;}
    void bind (Context *c) {}
    char *debug (char *buf) {sprintf (buf, "%g", value); return buf;}
};

class Scan;

class Named_Exp : public Expr  // Virtual class for names and funcs
{
protected:
    char *name;
public:
    virtual ~Named_Exp () {free (name);}
    Named_Exp (char *n) {name=n;}
    friend Named_Exp *Name_Or_Func_Expr (char *n, Scan *s);
};

Named_Exp *Name_Or_Func_Expr (char *n, Scan *s);

class Name_Expr : public Named_Exp
{
    NAME_FUNCP fptr;
    int key;
public:
    Name_Expr (char *n) : Named_Exp (n) {}
    void bind (Context *c);
    virtual double eval () 
	{if (!fptr) throw Unresolved_Name(); return (*fptr)(key);}
    char *debug (char *buf) {sprintf (buf,"[%s]",name);
                             return buf;}
};

struct function_node
{
    const char *name;
    int key;
    double (*fptr)(int*, double*, double*);
    function_node *next;
};

void function_setup (const char *name, int key, 
		       double (*fptr)(int*,double*,double*));

class Function_Expr : public Named_Exp
{
    Expr *arg;
    ~Function_Expr () {delete arg;}
    double (*fptr) (int*, double*, double*);
    static function_node *function_list;
    static function_node *it;
public:
    friend void function_setup (const char *name, int key, 
		       double (*fptr)(int*,double*,double*));
    Function_Expr (char *n, Expr *argument);
    double eval () {int i = 0; double a = arg->eval(); double d = 0.; 
                    return (*fptr)(&i, &a, &d);}
    static void list_start () {it = function_list;}
    static int list_ok () {return (0 != it);}
    static const char *list_next () {const char *t=it->name; it=it->next; 
                                     return t;}
    void bind (Context *c) {arg->bind(c);}
    char *debug (char *buf) {sprintf (buf, "%s(", name); 
                             arg->debug (&buf[strlen(buf)]);
			     strcat (buf, ")"); return buf;}


};

class Parenthesized_Expr : public Expr
{
public:
    Expr *expr;
    ~Parenthesized_Expr () {delete expr;}
    double eval () {return expr->eval();}
    void bind (Context *c) {expr->bind (c);}
    char *debug (char *buf) 
	{strcpy (buf,"("); expr->debug(&buf[1]); return strcat (buf, ")");}
};

class Binary_Expr : public Expr
{
public:
    Expr *first_expr;
    Expr *second_expr;
    ~Binary_Expr () {delete first_expr; delete second_expr;}
    virtual double eval () = 0;
    virtual Expr *insert (Expr *new_expr, Expr **top);
    void bind (Context *c) 
	{first_expr->bind (c); second_expr->bind (c);}
};

class Addition_Expr : public Binary_Expr
{
public:
    int precedence () {return Add_Precedence;}
    double eval () {return first_expr->eval() + second_expr->eval();}
    char *debug (char *buf) 
      {first_expr->debug (buf); strcat (buf, "+"); second_expr->debug 
       (&buf[strlen (buf)]); return buf;}
};

class Subtraction_Expr : public Binary_Expr
{
public:
    int precedence () {return Add_Precedence;}
    double eval () {return first_expr->eval() - second_expr->eval();}
    char *debug (char *buf) 
      {first_expr->debug (buf); strcat (buf, "-"); second_expr->debug 
       (&buf[strlen (buf)]); return buf;}
};

class Multiplication_Expr : public Binary_Expr
{
public:
    int precedence () {return Mult_Precedence;}
    double eval () {double d = first_expr->eval();
	return d ? d * second_expr->eval() : d;}
    char *debug (char *buf) 
      {first_expr->debug (buf); strcat (buf, "*"); second_expr->debug 
       (&buf[strlen (buf)]); return buf;}
};

class Division_Expr : public Binary_Expr
{
public:
    int precedence () {return Mult_Precedence;}
    double eval () {return first_expr->eval() / second_expr->eval();}
    char *debug (char *buf) 
      {first_expr->debug (buf); strcat (buf, "/"); second_expr->debug 
       (&buf[strlen (buf)]); return buf;}
};

class Power_Expr : public Binary_Expr
{
public:
    int precedence () {return Power_Precedence;}
    double eval () {return pow (first_expr->eval(), second_expr->eval());}
    char *debug (char *buf) 
      {first_expr->debug (buf); strcat (buf,"**"); second_expr->debug 
       (&buf[strlen (buf)]); return buf;}
};

class Eq_Expr : public Binary_Expr
{
public:
    int precedence () {return Eq_Precedence;}
    double eval () {return (first_expr->eval() == second_expr->eval());}
    char *debug (char *buf)
	{first_expr->debug (buf); strcat (buf,"=="); 
	second_expr->debug (&buf[strlen (buf)]); return buf;}
    bool boolean () {return true;}
};

class Ge_Expr : public Binary_Expr
{
public:
    int precedence () {return Eq_Precedence;}
    double eval () {return (first_expr->eval() >= second_expr->eval());}
    char *debug (char *buf)
	{first_expr->debug (buf); strcat (buf,">="); 
	second_expr->debug (&buf[strlen (buf)]); return buf;}
    bool boolean () {return true;}
};

class Le_Expr : public Binary_Expr
{
public:
    int precedence () {return Eq_Precedence;}
    double eval () {return (first_expr->eval() <= second_expr->eval());}
    char *debug (char *buf)
	{first_expr->debug (buf); strcat (buf,"<="); 
	second_expr->debug (&buf[strlen (buf)]); return buf;}
    bool boolean () {return true;}
};

class Ne_Expr : public Binary_Expr
{
public:
    int precedence () {return Eq_Precedence;}
    double eval () {return (first_expr->eval() != second_expr->eval());}
    char *debug (char *buf)
	{first_expr->debug (buf); strcat (buf,"!="); 
	second_expr->debug (&buf[strlen (buf)]); return buf;}
    bool boolean () {return true;}
};

class Gt_Expr : public Binary_Expr
{
public:
    int precedence () {return Eq_Precedence;}
    double eval () {return (first_expr->eval() > second_expr->eval());}
    char *debug (char *buf)
	{first_expr->debug (buf); strcat (buf,">"); 
	second_expr->debug (&buf[strlen (buf)]); return buf;}
    bool boolean () {return true;}
};

class Lt_Expr : public Binary_Expr
{
public:
    int precedence () {return Eq_Precedence;}
    double eval () {return (first_expr->eval() < second_expr->eval());}
    char *debug (char *buf)
	{first_expr->debug (buf); strcat (buf,"< "); 
	second_expr->debug (&buf[strlen (buf)]); return buf;}
    bool boolean () {return true;}
};


#endif



