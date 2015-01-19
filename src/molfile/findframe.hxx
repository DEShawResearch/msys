#include <cmath>
#include <sys/types.h>

namespace desres { namespace molfile { namespace findframe {

    /* The Oracle in all these functions implements operator[] taking a
     * ssize_t in the half-open range [0,N) and returning a time. */

    template <typename Oracle>
    void lookup_time( ssize_t N, const Oracle& times, double T, 
                         ssize_t& left, ssize_t& right ) {

       ssize_t mid;

       // -----------------------------------------------
       // Time not found if we have no frames!
       // -----------------------------------------------
       left = 0;
       right = N-1;
       if (N <= 0) return;

       // -----------------------------------------------
       // Knuth's binary search
       // Notice that mid = (left+right)/2 suffers from
       // overflow problems, so the less intuitive
       // mid = left + ((right-left)/2)
       // is used instead.
       // -----------------------------------------------
       while(left <= right) {
           mid = left + ((right-left)/2);
           double mid_time = times[mid];
           if (T > mid_time) {
               left = mid + 1;
           } else if (T < mid_time) {
               right = mid - 1;
           } else {
               /* Exact match */
               left = right = mid;
               break;
           }
       }

       // -----------------------------------------------
       // CAUTION: at this point, the meanings of 
       // left and right are switched (i.e. left >= right,
       // see the while() loop above if you don't believe me!
       // -----------------------------------------------
       ssize_t swap = left;
       left = right;
       right = swap;
    }   

    template <typename Oracle>
    ssize_t at_time_near( ssize_t N, const Oracle& times, double T ) {
        ssize_t left, right;
        lookup_time(N,times,T,left,right);
        if (right<0) return -1;     /* no frames */
        if (left<0) return 0;       /* time less than first time */
        if (right==N) return left;  /* right greater than last time */
        return (fabs(T-times[right]) < fabs(T-times[left])) ? right : left;
    }
    
    template <typename Oracle>
    ssize_t at_time_le(   ssize_t N, const Oracle& times, double T ) {
        ssize_t left, right;
        lookup_time(N,times,T,left,right);
        if (right<0) return -1;
        if (left<0) return -1;
        if (right>N) return -1;
        if (times[left]>T) return -1;
        return left;
    }
    
    template <typename Oracle>
    ssize_t at_time_lt(   ssize_t N, const Oracle& times, double T ) {
        ssize_t left, right;
        lookup_time(N,times,T,left,right);
        if (right<0) return -1;
        if (left<=0) return -1;
        if (left==right) {
            if (times[left-1]>=T) return -1;
            return left-1;
        }
        if (times[left]>=T) return -1;
        return left;
    }
    
    template <typename Oracle>
    ssize_t at_time_ge(   ssize_t N, const Oracle& times, double T ) {
        ssize_t left, right;
        lookup_time(N,times,T,left,right);
        if (right<0) return -1;
        if (right>=N) return -1;
        if (times[right]<T) return -1;
        return right;
    }
    
    template <typename Oracle>
    ssize_t at_time_gt(   ssize_t N, const Oracle& times, double T ) {
        ssize_t left, right;
        lookup_time(N,times,T,left,right);
        if (right<0) return -1;
        if (right>=N) return -1;
        if (left==right) {
            if (right+1==N) {
                if (times[left]>T) return left; /* single frame case */
                return -1;
            }
            if (times[right+1]<=T) return -1;
            return right+1;
        }
        if (times[right]<=T) return -1;
        return right;
    }

}}}
