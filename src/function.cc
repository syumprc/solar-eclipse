/*
 * function.cc sets up default intrinsic functions for expression.cc
 * Written by Charles Peterson beginning on October 21, 1997
 * Copyright (c) 1997 Southwest Foundation for Biomedical Research
 */

#include <math.h>
#include "token.h"
#include "expression.h"

// Note: gamma() has been changed to tgamma() for C++0x compatibility
// results seem to be the identical, though the original gamma() was
// not part of a standard and could have actually implemented log gamma
// on some systems


double dprint (int *key, double *argument, double *argument2)
{
    char buf[1024];
    printf ("%g >", *argument);
    fgets (buf, 1024, stdin);
    return *argument;
}

double dsqrt (int *key, double *argument, double *argument2)
{
    return sqrt (*argument);
}

double dlog (int *key, double *argument, double *argument2)
{
    return log (*argument);
}

double dlog10 (int *key, double *argument, double *argument2)
{
    return log10 (*argument);
}

double dabs (int *key, double *argument, double *argument2)
{
    return fabs (*argument);
}

double derf (int *key, double *argument, double *argument2)
{
    return erf (*argument);
}

double derfc (int *key, double *argument, double *argument2)
{
    return erfc (*argument);
}

double dgamma (int *key, double *argument, double *argument2)
{
    return tgamma (*argument);
}

double dlgamma (int *key, double *argument, double *argument2)
{
    return lgamma (*argument);
}

double dj0 (int *key, double *argument, double *argument2)
{
    return j0 (*argument);
}

double dj1 (int *key, double *argument, double *argument2)
{
    return j1 (*argument);
}

double dy0 (int *key, double *argument, double *argument2)
{
    return y0 (*argument);
}

double dy1 (int *key, double *argument, double *argument2)
{
    return y1 (*argument);
}

double dceil (int *key, double *argument, double *argument2)
{
    return ceil (*argument);
}

double dfloor (int *key, double *argument, double *argument2)
{
    return floor (*argument);
}

double drint (int *key, double *argument, double *argument2)
{
    return rint (*argument);
}

double dlog1p (int *key, double *argument, double *argument2)
{
    return log1p (*argument);
}

double dexpm1 (int *key, double *argument, double *argument2)
{
    return expm1 (*argument);
}

double dlogb (int *key, double *argument, double *argument2)
{
    return logb (*argument);
}

double dexp (int *key, double *argument, double *argument2)
{
    return exp (*argument);
}

double dcbrt (int *key, double *argument, double *argument2)
{
    return cbrt (*argument);
}

double dsin (int *key, double *argument, double *argument2)
{
    return sin (*argument);
}

double dcos (int *key, double *argument, double *argument2)
{
    return cos (*argument);
}

double dtan (int *key, double *argument, double *argument2)
{
    return tan (*argument);
}

double dasin (int *key, double *argument, double *argument2)
{
    return asin (*argument);
}

double dacos (int *key, double *argument, double *argument2)
{
    return acos (*argument);
}

double datan (int *key, double *argument, double *argument2)
{
    return atan (*argument);
}

double dsinh (int *key, double *argument, double *argument2)
{
    return sinh (*argument);
}

double dcosh (int *key, double *argument, double *argument2)
{
    return cosh (*argument);
}

double dtanh (int *key, double *argument, double *argument2)
{
    return tanh (*argument);
}

void setup_functions ()
{
    function_setup ("erfc", 0, derfc);
    function_setup ("erf", 0, derf);
    function_setup ("lgamma", 0, dlgamma);
    function_setup ("gamma", 0, dgamma);
    function_setup ("j1", 0, dj1);
    function_setup ("j0", 0, dj0);
    function_setup ("y1", 0, dy1);
    function_setup ("y0", 0, dy0);
    function_setup ("rint", 0, drint);
    function_setup ("floor", 0, dfloor);
    function_setup ("ceil", 0, dceil);
    function_setup ("tanh", 0, dtanh);
    function_setup ("cosh", 0, dcosh);
    function_setup ("sinh", 0, dsinh);
    function_setup ("atan", 0, datan);
    function_setup ("acos", 0, dacos);
    function_setup ("asin", 0, dasin);
    function_setup ("tan", 0, dtan);
    function_setup ("cos", 0, dcos);
    function_setup ("sin", 0, dsin);
    function_setup ("expm1", 0, dexpm1);
    function_setup ("exp", 0, dexp);
    function_setup ("logb", 0, dlogb);
    function_setup ("log1p", 0, dlog1p);
    function_setup ("log10", 0, dlog10);
    function_setup ("log", 0, dlog);
    function_setup ("cbrt", 0, dcbrt);
    function_setup ("sqrt", 0, dsqrt);
    function_setup ("abs", 0, dabs);
    function_setup ("print", 0, dprint);
}


