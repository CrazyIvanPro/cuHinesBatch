/**
 *
 * 	@file utils.h

 *
 * 	@brief utils funtions
 *
 * 	@author Ivan Martinez-Perez ivan.martinez@bsc.es
 * 	@author Pedro Valero-Lara   pedro.valero@bsc.es
 *
 **/
#include <random>

#if defined( _WIN32 ) || defined( _WIN64 )
#  include <time.h>
#  include <sys/timeb.h>
#  if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#    define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#  else
#    define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#  endif
#else
#  include <sys/time.h>
#  include <time.h>
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- Emulate gettimeofday on Windows.
*/
#if defined( _WIN32 ) || defined( _WIN64 )
#ifndef _TIMEZONE_DEFINED
#define _TIMEZONE_DEFINED

struct timezone
{
    int  tz_minuteswest; /* minutes W of Greenwich */
    int  tz_dsttime;     /* type of dst correction */
};
#endif

int gettimeofday(struct timeval* tv, struct timezone* tz)
{
    FILETIME         ft;
    unsigned __int64 tmpres = 0;
    static int       tzflag = 0;

    if (NULL != tv) {
        GetSystemTimeAsFileTime(&ft);
        tmpres |=  ft.dwHighDateTime;
        tmpres <<= 32;
        tmpres |=  ft.dwLowDateTime;

        /*converting file time to unix epoch*/
        tmpres /= 10;  /*convert into microseconds*/
        tmpres -= DELTA_EPOCH_IN_MICROSECS;

        tv->tv_sec  = (long)(tmpres / 1000000UL);
        tv->tv_usec = (long)(tmpres % 1000000UL);
    }
    if (NULL != tz) {
        if (!tzflag) {
            _tzset();
            tzflag = 1;
        }
        tz->tz_minuteswest = _timezone / 60;
        tz->tz_dsttime     = _daylight;
    }
    return 0;
}
#endif


/**
 *  Return time in seconds since arbitrary time in the past.
 *  Used for elapsed wall clock time computation.
 **/
double time_wtime( void )
{
    struct timeval t;
    gettimeofday( &t, NULL );
    return t.tv_sec + t.tv_usec*1e-6;
}

/**
 * Version callable from Fortran stores seconds in the parameter time.
 **/
void timef_wtime(double *time)
{
    *time = time_wtime();
}


void check(){
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) 
        printf("Error: %s\n", cudaGetErrorString(err));
}


double doubleRand(double fMin, double fMax){

    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);

}

extern int integerRand(int iMin, int iMax){
    return rand()%(iMax-iMin + 1) + iMin;
}

void calcError(double *src,double *dst, int size){

    double error=0;

    for (int i = 0; i < size; ++i)
    {
		//printf("src[%i]=%.2f-dts[%i]=%.2f\n",i,src[i],i,dst[i]);

        if (error < abs(abs(src[i]) - abs(dst[i])))
        {
            error = abs(abs(src[i]) - abs(dst[i]));
        }

        
    }
    printf("    error: %e\n",error);
}
