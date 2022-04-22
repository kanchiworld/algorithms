/*
This is a demo code for AKS primality testing algorithm

AKS is NOT the fastest primality testing algorithm. For practical usage better algorithms are avaialble. The value of AKS is
primarily theoretical. It was the first algo to prove that primality test lies in the polynomial complexity class.

This is also a very novice implementation of AKS. long integer type is used for input. For an industrial implementation one will
need some big integer class. On a 64-bit machine the code conks out for integers larger than 10 digits. However, the novice
implementation also makes the code easy to understand.

compile as 
g++ -std=c++11 aks.cpp -o aks

runas as
aks <input num>
for e.g
aks 1973

*/

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>

// a^e mod m
unsigned long modexp(unsigned long a, unsigned long e, unsigned long m);

//Euler's totient function
unsigned long phi(unsigned long n);

//gcd
long gcd(long n1, long n2);

// poly1 * poly2
template <typename T>
void PolyMultiply(const std::vector<T>& poly1, const std::vector<T>& poly2, std::vector<T>& output);

// basePoly^e mod moduloPoly mod m
template <typename T>
void PolyModExp(const std::vector<T>& basePoly, unsigned long e, const std::vector<T>& moduloPoly, unsigned long m, std::vector<T>& remainPoly);

// numer/denom
template <typename T>
void PolyDivide(const std::vector<T>& numer, const std::vector<T>& denom, std::vector<T>& quot, std::vector<T>& remain);

// poly1 - poly2
template <typename T>
void PolySubtract(const std::vector<T>& poly1, const std::vector<T>& poly2, std::vector<T>& output);

template <typename T>
void PrintV(const std::string& name, const std::vector<T>& A);


int main(int argc, char *argv[])
{
    if ( argc != 2 )
    {
        std::cout << "Usage: aks <positive integer>" << std::endl;
        exit(-1);
    }

    unsigned long n = atol(argv[1]); 
//  std::cout << "n = " << n << std::endl ;

//  Find the smallest r such that Or(n) > (log2 n)2.
//  maxk=⌊(log2 n)2⌋;
    unsigned long ceillog2n = 0, tempn = n;
    while ( tempn ) { tempn >>= 1; ceillog2n++;} ;
//  ceillog2n will give correct result for all postive ints except those of the form 2^k, in which case it will overestimate them by a count of 1. not a big deal.
//  std::cout << "ceillog2n = " << ceillog2n << std::endl ;
    unsigned long maxk = (ceillog2n) * (ceillog2n) ;

//  maxr=Max[3, ⌈(Log2 n)5⌉]; (*maxr really isn't needed*)
    unsigned long maxr = std::max(3UL, ceillog2n*ceillog2n*ceillog2n*ceillog2n*ceillog2n); 
    bool nextR = true;
    unsigned long r = 2;
    for (; nextR && r < maxr; r++)
    {
        nextR = false;
        for (unsigned long k = 1; !nextR && k <= maxk; k++)
        {
            unsigned long mexp = modexp(n, k, r);
//          std::cout << "[mexp = modexp(n, k, r)] " << mexp << " = modexp(" << n << ", " << k << ", " << r << ")" << std::endl;
            nextR = (mexp == 1 || mexp == 0);
        }
    };
    r--; //(*the loop over increments by one*)
//    std::cout << "r = " << r << std::endl; 

//  If 1 < gcd(a,n) < n for some a ≤ r, output composite.
    for (long a = r; a > 1; a--)
    {
        long g = gcd(a,n);
        if (g > 1 && g < n)
        {
            std::cout <<   "COMPOSITE"  << std::endl;
            return 0;
        }
    }
//    gcd={GCD(29,31)=1, GCD(28,31)=1, ..., GCD(2,31)=1} ≯ 1
     
//  If n ≤ r, output prime.
    if (n <= r)
    {
            std::cout <<   "PRIME"  << std::endl;
            return 0;
// (* this step may be omitted if n > 5690034 *)
//    31 > 29
    }
     
    unsigned long fi = phi(r);
//  std::cout << "\nfi = " << fi << std::endl; 

    unsigned long maxa = (unsigned long)floor(sqrt(fi)*ceillog2n);
//  std::cout << "\nmaxa = " << maxa << std::endl; 

    std::vector<long> moduloPoly(r+1, 0);
    moduloPoly[r] = 1;
    moduloPoly[0] = -1;

    std::vector<long> lhsRemainPoly;
    std::vector<long> rhsRemainPoly(n%r+1, 0);
    std::vector<long> remainPoly;

    bool isPrime = true;    
    for (long a = 1; a <= maxa; a++)
    {
        const std::vector<long> basePoly{a, 1};
        PolyModExp(basePoly, n, moduloPoly, n, lhsRemainPoly);
      
        rhsRemainPoly[n%r] = 1; rhsRemainPoly[0] = a;

        PolySubtract(lhsRemainPoly, rhsRemainPoly, remainPoly);

        if ( remainPoly.size() > 0 )
        {
            isPrime = false;
            break;
        }
    }
    std::cout << ( isPrime ? "PRIME" :  "COMPOSITE" ) << std::endl;

    return 0;
}

//modular exponentiation function. returns a^e mod m
unsigned long modexp(unsigned long a, unsigned long e, unsigned long m)
{
    unsigned long retval = 1;
    while (e)
    {
        if ( e & 0x01UL ) retval = retval*a % m;
        e >>= 1;
        if ( e ) a = (a*a)%m ;
    }
    return retval;
}

//gcd
long gcd(long n1, long n2)
{
   long r = 1;
    while( (r = n1%n2) )
   {
       n1 = n2;
       n2 = r;
   }
   return n2;
}



//Euler's totient function
unsigned long phi(unsigned long n)
{
    unsigned long result = n;
    for (unsigned long i = 2; i * i <= n; i++)
    {
        if (n % i == 0)
        {
            while (n % i == 0) n /= i;
            result -= result / i;  //this is n(1-1/p1) = n - n/p1
        }
    }
    if (n > 1) result -= result / n;
    return result;
}


 
 
// prints all members of the vector
template <typename T>
void PrintV(const std::string& name, const std::vector<T>& A)
{
    std::cout << name << "(size=" << A.size() << "), val[0..n] = [ ";
    for (auto const &i: A) { std::cout << i << " "; }
    std::cout << "]\n";
}


 
//Calculates polynomial division numer/denom
template <typename T>
void PolyDivide(const std::vector<T>& numer, const std::vector<T>& denom, std::vector<T>& quot, std::vector<T>& remain)
{
    // vectors - N / D == q && N % D == r
//  make a copy of numer into remain 
    remain = numer;

    if ( denom.size() > numer.size() ) //Nothing to be done
        return;

    quot.resize(numer.size()-denom.size()+1);
    int dN = numer.size()-1;
    int dD = denom.size()-1;

    while( dN >= dD  )
    {
        // calculating one element of q
        quot[dN-dD] = remain[dN]/denom[dD];
 
        for( int i = 0 ; i <= dD ; i++ )
            remain[dN-i] -= denom[dD-i] * quot[dN-dD];
 
        dN--;
    }
    while (! remain.back()) remain.pop_back();
}

//Calculates poly1-poly2
template <typename T>
void PolySubtract(const std::vector<T>& poly1, const std::vector<T>& poly2, std::vector<T>& output)
{
    output.resize( std::max(poly1.size(), poly2.size()));
    unsigned long commonlength = std::min(poly1.size(), poly2.size());

    for ( int i = 0; i <= commonlength-1; i++)
            output[i] = poly1[i] - poly2[i];

    int sign = poly1.size() >= poly2.size() ? 1 : -1 ;

    for ( int i = commonlength; i <= output.size()-1; i++)
        output[i] = sign == 1 ? poly1[i] : -1 * poly2[i];

    while (! output.back()) output.pop_back();
}

 
//Calculates poly1*poly2
template <typename T>
void PolyMultiply(const std::vector<T>& poly1, const std::vector<T>& poly2, std::vector<T>& output)
{
    output.resize(poly1.size()+poly2.size()-1, 0);

    for ( int i = 0; i <= poly1.size()-1; i++)
        for ( int j = 0; j <= poly2.size()-1; j++)
            output[j+i] += poly2[j]*poly1[i];
}

 
//polynomial modular exponentiation function. returns ( basePoly^e mod moduloPoly) mod m
template <typename T>
void PolyModExp(const std::vector<T>& basePoly, unsigned long e, const std::vector<T>& moduloPoly, unsigned long m, std::vector<T>& remainPoly)
{
    std::vector<T> a = basePoly;

    remainPoly.resize(1);
    remainPoly[0] = 1;

    while (e)
    {
        std::vector<T> tempmult;
        std::vector<T> q; //throw away val
        if ( e & 0x01UL )
        {
            PolyMultiply(remainPoly, a, tempmult) ;
            std::for_each(tempmult.begin(), tempmult.end(), [m](T& n){ n %= m; });
            PolyDivide(tempmult, moduloPoly, q, remainPoly);
        }
        e >>= 1;
        if ( e )
        {
        std::vector<T> asq;
        PolyMultiply(a, a, asq);
        std::for_each(asq.begin(), asq.end(), [m](T& n){ n %= m; });
        PolyDivide(asq, moduloPoly, q, a);
        }
    }
    std::for_each(remainPoly.begin(), remainPoly.end(), [m](T& n){ n %= m; });
    while (! remainPoly.back()) remainPoly.pop_back();
    //PrintV("remainPoly", remainPoly);
}

