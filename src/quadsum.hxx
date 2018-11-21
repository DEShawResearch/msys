#ifndef desres_msys_quadsum_hxx
#define desres_msys_quadsum_hxx

namespace desres { namespace msys {

    // compensated summation (sumXBLAS from "ACCURATE SUM AND DOT PRODUCT")
    class Quadsum {
    public:
        Quadsum(): total(0), err(0.0){};
        Quadsum(double t): total(t),err(0.0){};
        Quadsum(double t, double e): total(t),err(e){};
        Quadsum(const Quadsum &qs): total(qs.total),err(qs.err){};
    
        Quadsum & operator+=(const double val){
            double t1,t2;
            t1=twosum(total,val,t2);
            t2+=err;
            total=fasttwosum(t1,t2,err);
            return *this;
        }
        Quadsum & operator+=(const Quadsum val){
            double t1,t2;
            t1=twosum(total,val.total,t2);
            t2+=err+val.err;
            total=fasttwosum(t1,t2,err);
            return *this;
        }
        Quadsum operator+(const double val) const{
            Quadsum tmp(total,err);
            tmp+=val;
            return tmp;
        }
        Quadsum operator+(const Quadsum val) const{
            Quadsum tmp(total,err);
            tmp+=val;
            return tmp;
        }
    
        Quadsum & operator=(const double val){
            total=val;
            err=0;
            return *this;
        }
        Quadsum & operator=(const Quadsum &val){
            total=val.total;
            err=val.err;
            return *this;
        }
        double result(){return total;}
        void reset(){total=0.0;err=0.0;}
    private:
        double total, err;
        /* "ACCURATE SUM AND DOT PRODUCT", SIAM Journal on 
           Scientific Computing 26(6):1955-1988, 2005 */
        double twosum(const double a, const double b, double &y){
            double xt,zt;
            xt=a+b;
            zt=xt-a;
            y=(a-(xt-zt)+(b-zt));
            return xt;
        }
        double fasttwosum(const double a, const double b, double &y){
            double xt=a+b;
            y=(a-xt)+b;
            return xt;
        }
    };

}}

#endif

