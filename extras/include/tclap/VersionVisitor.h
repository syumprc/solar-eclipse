// -*- Mode: c++; c-basic-offset: 4; tab-width: 4; -*-

/****************************************************************************** 
 * 
 *  file:  VersionVisitor.h
 * 
 *  
 *****************************************************************************/ 


#ifndef TCLAP_VERSION_VISITOR_H
#define TCLAP_VERSION_VISITOR_H

#include <tclap/CmdLineInterface.h>
#include <tclap/CmdLineOutput.h>
#include <tclap/Visitor.h>

namespace TCLAP {

/**
 * A Vistor that will call the version method of the given CmdLineOutput
 * for the specified CmdLine object and then exit.
 */
class VersionVisitor: public Visitor
{
	protected:

		/**
		 * The CmdLine of interest.
		 */
		CmdLineInterface* _cmd;

		/**
		 * The output object. 
		 */
		CmdLineOutput** _out;

	public:

		/**
		 * Constructor.
		 * \param cmd - The CmdLine the output is generated for. 
		 * \param out - The type of output. 
		 */
		VersionVisitor( CmdLineInterface* cmd, CmdLineOutput** out ) 
				: Visitor(), _cmd( cmd ), _out( out ) { }

		/**
		 * Calls the version method of the output object using the
		 * specified CmdLine.
		 */
		void visit() { 
		    (*_out)->version(*_cmd); 
		    throw ExitException(0); 
		}

};

}

#endif
