		#ifndef NJAKSDTHASKJERAXJGFBZJDLAGZ
#define NJAKSDTHASKJERAXJGFBZJDLAGZ

// At most one of the flags
//    PM_TIMER_MSVC
//    PM_TIMER_CLOCK_GETTIME
//    PM_TIMER_GETRUSAGE
//    PM_TIMER_EXTERNAL
//    PM_TIMER_NONE
// can be defined externally. If PM_TIMER_EXTERNAL is defined,
// then there must exist a definition of function "double get_time()" elsewhere.

#if defined (PM_TIMER_MSVC) || defined (PM_TIMER_CLOCK_GETTIME) || defined (PM_TIMER_GETRUSAGE) || defined (PM_TIMER_EXTERNAL) || defined (PM_TIMER_NONE)
#else
	// default option
	#ifdef _MSC_VER
		#define PM_TIMER_MSVC
	#else
		#define PM_TIMER_CLOCK_GETTIME
	#endif
#endif

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

#ifdef PM_TIMER_MSVC

	#include <windows.h>

	inline double get_time()
	{
		LARGE_INTEGER t, frequency;
		QueryPerformanceCounter(&t);
		QueryPerformanceFrequency(&frequency);
		return (double)t.QuadPart/(double)frequency.QuadPart;
	}

#endif

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

#ifdef PM_TIMER_CLOCK_GETTIME

//	#include <time.h>
//
//	inline double get_time()
//	{
//		struct timespec t;
//		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
//		return (double)t.tv_nsec*1.00E-9 + (double)t.tv_sec;
//	}


#include <time.h>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif


inline double get_time()
{
    struct timespec t;
    
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    t.tv_sec = mts.tv_sec;
    t.tv_nsec = mts.tv_nsec;
#else
    clock_gettime(CLOCK_REALTIME, &t);
#endif
    
    return (double)t.tv_nsec*1.00E-9 + (double)t.tv_sec;
    
}



#endif

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

#ifdef PM_TIMER_GETRUSAGE

	#include <sys/resource.h>

	inline double get_time()
	{
		struct rusage t;
		getrusage (RUSAGE_SELF, &t);
		return (double)t.ru_utime.tv_usec*1.00E-6 + (double)t.ru_utime.tv_sec;
	}

#endif

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

#ifdef PM_TIMER_EXTERNAL

	extern double get_time();

#endif

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

#ifdef PM_TIMER_NONE

	inline double get_time() { return 0; }

#endif

#endif



