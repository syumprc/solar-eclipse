#include "solar.h"

extern "C" int Solar_Init (Tcl_Interp *interp);

extern "C" int create_interp(ClientData clientData, Tcl_Interp *interp,
                                         int argc,const char *argv[]){
											 
	Tcl_FindExecutable( NULL ); 										 
	Tcl_Interp * interp_2 = Tcl_CreateInterp ();
	Tcl_Init(interp_2);
	
	Solar_Init(interp_2);
	
	Tcl_Eval(interp_2, "trait trt_2_0001");
	
	
	
	return TCL_OK;
}
