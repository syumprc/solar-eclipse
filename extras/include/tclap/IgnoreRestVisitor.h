
/****************************************************************************** 
 * 
 *  file:  IgnoreRestVisitor.h
 *  
 *  
 *****************************************************************************/ 


#ifndef TCLAP_IGNORE_REST_VISITOR_H
#define TCLAP_IGNORE_REST_VISITOR_H

#include <tclap/Visitor.h>
#include <tclap/Arg.h>

namespace TCLAP {

/**
 * A Vistor that tells the CmdLine to begin ignoring arguments after
 * this one is parsed.
 */
class IgnoreRestVisitor: public Visitor
{
	public:

		/**
		 * Constructor.
		 */
		IgnoreRestVisitor() : Visitor() {}

		/**
		 * Sets Arg::_ignoreRest.
		 */
		void visit() { Arg::beginIgnoring();  }
};

}

#endif
