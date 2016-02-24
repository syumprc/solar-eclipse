

/****************************************************************************** 
 * 
 *  file:  CmdLineOutput.h
 * 
 *  
 *****************************************************************************/ 

#ifndef TCLAP_CMDLINEOUTPUT_H
#define TCLAP_CMDLINEOUTPUT_H

#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace TCLAP {

class CmdLineInterface;
class ArgException;

/**
 * The interface that any output object must implement.
 */
class CmdLineOutput 
{

	public:

		/**
		 * Virtual destructor.
		 */
		virtual ~CmdLineOutput() {}

		/**
		 * Generates some sort of output for the USAGE. 
		 * \param c - The CmdLine object the output is generated for. 
		 */
		virtual void usage(CmdLineInterface& c)=0;

		/**
		 * Generates some sort of output for the version. 
		 * \param c - The CmdLine object the output is generated for. 
		 */
		virtual void version(CmdLineInterface& c)=0;

		/**
		 * Generates some sort of output for a failure. 
		 * \param c - The CmdLine object the output is generated for. 
		 * \param e - The ArgException that caused the failure. 
		 */
		virtual void failure( CmdLineInterface& c, 
						      ArgException& e )=0;

};

} //namespace TCLAP
#endif 
