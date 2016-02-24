
/****************************************************************************** 
 * 
 *  file:  HelpVisitor.h
 * 
 *  
 *****************************************************************************/ 

#ifndef TCLAP_HELP_VISITOR_H
#define TCLAP_HELP_VISITOR_H

#include <tclap/CmdLineInterface.h>
#include <tclap/CmdLineOutput.h>
#include <tclap/Visitor.h>

namespace TCLAP {

/**
 * A Visitor object that calls the usage method of the given CmdLineOutput
 * object for the specified CmdLine object.
 */
class HelpVisitor: public Visitor
{
	protected:

		/**
		 * The CmdLine the output will be generated for. 
		 */
		CmdLineInterface* _cmd;

		/**
		 * The output object. 
		 */
		CmdLineOutput** _out;

	public:

		/**
		 * Constructor.
		 * \param cmd - The CmdLine the output will be generated for.
		 * \param out - The type of output. 
		 */
		HelpVisitor(CmdLineInterface* cmd, CmdLineOutput** out) 
				: Visitor(), _cmd( cmd ), _out( out ) { }

		/**
		 * Calls the usage method of the CmdLineOutput for the 
		 * specified CmdLine.
		 */
		void visit() { (*_out)->usage(*_cmd); throw ExitException(0); }
		
};

}

#endif
