
/****************************************************************************** 
 * 
 *  file:  Visitor.h
 * 
 *  
 *****************************************************************************/ 


#ifndef TCLAP_VISITOR_H
#define TCLAP_VISITOR_H

namespace TCLAP {

/**
 * A base class that defines the interface for visitors.
 */
class Visitor
{
	public:

		/**
		 * Constructor. Does nothing.
		 */
		Visitor() { }

		/**
		 * Destructor. Does nothing.
		 */
		virtual ~Visitor() { }

		/**
		 * Does nothing. Should be overridden by child.
		 */
		virtual void visit() { }
};

}

#endif
