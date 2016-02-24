/*
*
*
*/



#ifndef STOPWATCH_H
#define STOPWATCH_H

// for clock() and CLOCKS_PER_SEC
#include <time.h>


namespace TNT
{

inline static double seconds(void)
{
    const double secs_per_tick = 1.0 / CLOCKS_PER_SEC;
    return ( (double) clock() ) * secs_per_tick;
}

class Stopwatch {
    private:
        int running_;
        double start_time_;
        double total_;

    public:
        inline Stopwatch();
        inline void start();
        inline double stop();
		inline double read();
		inline void resume();
		inline int running();
};

inline Stopwatch::Stopwatch() : running_(0), start_time_(0.0), total_(0.0) {}

void Stopwatch::start() 
{
	running_ = 1;
	total_ = 0.0;
	start_time_ = seconds();
}

double Stopwatch::stop()  
{
	if (running_) 
	{
         total_ += (seconds() - start_time_); 
         running_ = 0;
    }
    return total_; 
}

inline void Stopwatch::resume()
{
	if (!running_)
	{
		start_time_ = seconds();
		running_ = 1;
	}
}
		

inline double Stopwatch::read()   
{
	if (running_)
	{
		stop();
		resume();
	}
	return total_;
}


} /* TNT namespace */
#endif
    

            
