

/****************************************************************************** 
 * 
 *  file:  OptionalUnlabeledTracker.h
 * 
 *  
 *****************************************************************************/ 


#ifndef TCLAP_OPTIONAL_UNLABELED_TRACKER_H
#define TCLAP_OPTIONAL_UNLABELED_TRACKER_H

#include <string>

namespace TCLAP {

class OptionalUnlabeledTracker
{

	public:

		static void check( bool req, const std::string& argName );

		static void gotOptional() { alreadyOptionalRef() = true; }

		static bool& alreadyOptional() { return alreadyOptionalRef(); } 

	private:

		static bool& alreadyOptionalRef() { static bool ct = false; return ct; }
};


inline void OptionalUnlabeledTracker::check( bool req, const std::string& argName )
{
    if ( OptionalUnlabeledTracker::alreadyOptional() )
        throw( SpecificationException(
	"You can't specify ANY Unlabeled Arg following an optional Unlabeled Arg",
	                argName ) );

    if ( !req )
        OptionalUnlabeledTracker::gotOptional();
}


} // namespace TCLAP

#endif
