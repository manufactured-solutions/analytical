#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "nsctpl_fwd.hpp"
#include "nsctpl.hpp"

// Explicitly instantiate the primitive solution to ensure compilation
template class nsctpl::primitive_solution<float>;
template class nsctpl::primitive_solution<double>;
template class nsctpl::primitive_solution<long double>;

// Explicitly instantiate the manufactured forcing to ensure compilation
template class nsctpl::manufactured_solution<float>;
template class nsctpl::manufactured_solution<double>;
template class nsctpl::manufactured_solution<long double>;

// Used to output each solution parameter on std::cout
template<typename Scalar>
void print_it(const std::string &name, Scalar value) {
    std::cout << '\t' << name << " = " << value << std::endl;
}

// Data employed within prime_it below.  Taken from `primes 1 2000`
static const unsigned int primes[303] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151,
    157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233,
    239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317,
    331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419,
    421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
    509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607,
    613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701,
    709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811,
    821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911,
    919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013,
    1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091,
    1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181,
    1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277,
    1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361,
    1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451,
    1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531,
    1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609,
    1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699,
    1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789,
    1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
    1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997,
    1999
};

// Stateful functor used to set parameter values to the prime number sequence
template<typename Scalar>
struct prime_it {
    std::size_t i;
    prime_it() : i(0) {}
    void operator()(const std::string&, Scalar& value) {
        assert(i < sizeof(primes)/sizeof(primes[0]));
        value = primes[i++];
    }

    // Non-copyable as copying ruins the stateful behavior
private:
    prime_it (const prime_it &);
    prime_it& operator=(const prime_it &);
};

// Helper macro to dump results to std::cout
#define OUTPUT(expr) std::cout << '\t' << #expr << " = " << (expr) << std::endl;

int main(int argc, char *argv[])
{
    nsctpl::manufactured_solution<long double> ms;
    {
        prime_it<long double> p;
        ms.foreach_parameter(p);
    }

    std::cout << "Parameters after initialization:" << std::endl;
    ms.foreach_parameter(print_it<long double>);

    std::cout.precision(std::numeric_limits<long double>::digits10);

    std::cout << "Analytic solution evaluations:" << std::endl;
    OUTPUT(ms.soln_rho    (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_rho._t (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_rho._x (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_rho._xx(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_rho._xy(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_rho._xz(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_rho._y (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_rho._yy(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_rho._yz(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_rho._z (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_rho._zz(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_u    (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_u._t (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_u._x (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_u._xx(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_u._xy(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_u._xz(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_u._y (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_u._yy(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_u._yz(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_u._z (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_u._zz(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_v    (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_v._t (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_v._x (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_v._xx(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_v._xy(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_v._xz(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_v._y (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_v._yy(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_v._yz(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_v._z (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_v._zz(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_w    (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_w._t (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_w._x (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_w._xx(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_w._xy(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_w._xz(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_w._y (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_w._yy(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_w._yz(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_w._z (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_w._zz(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_T    (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_T._t (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_T._x (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_T._xx(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_T._xy(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_T._xz(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_T._y (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_T._yy(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_T._yz(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_T._z (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.soln_T._zz(0.5, 0.6, 0.7, 0.8));

    std::cout << "Analytically determined quantities" << std::endl;
    OUTPUT(ms.eval_exact_rho(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.eval_exact_u  (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.eval_exact_v  (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.eval_exact_w  (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.eval_exact_T  (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.eval_g_rho(0.5, 0.6, 0.7, 0.8, 1));
    OUTPUT(ms.eval_g_rho(0.5, 0.6, 0.7, 0.8, 2));
    OUTPUT(ms.eval_g_rho(0.5, 0.6, 0.7, 0.8, 3));
    OUTPUT(ms.eval_g_u  (0.5, 0.6, 0.7, 0.8, 1));
    OUTPUT(ms.eval_g_u  (0.5, 0.6, 0.7, 0.8, 2));
    OUTPUT(ms.eval_g_u  (0.5, 0.6, 0.7, 0.8, 3));
    OUTPUT(ms.eval_g_v  (0.5, 0.6, 0.7, 0.8, 1));
    OUTPUT(ms.eval_g_v  (0.5, 0.6, 0.7, 0.8, 2));
    OUTPUT(ms.eval_g_v  (0.5, 0.6, 0.7, 0.8, 3));
    OUTPUT(ms.eval_g_w  (0.5, 0.6, 0.7, 0.8, 1));
    OUTPUT(ms.eval_g_w  (0.5, 0.6, 0.7, 0.8, 2));
    OUTPUT(ms.eval_g_w  (0.5, 0.6, 0.7, 0.8, 3));
    OUTPUT(ms.eval_g_T  (0.5, 0.6, 0.7, 0.8, 1));
    OUTPUT(ms.eval_g_T  (0.5, 0.6, 0.7, 0.8, 2));
    OUTPUT(ms.eval_g_T  (0.5, 0.6, 0.7, 0.8, 3));

    std::cout << "Quantities built from the analytical solutions" << std::endl;
    OUTPUT(ms.eval_exact_e (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.eval_exact_p (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.eval_exact_mu(0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.eval_g_e     (0.5, 0.6, 0.7, 0.8, 1));
    OUTPUT(ms.eval_g_e     (0.5, 0.6, 0.7, 0.8, 2));
    OUTPUT(ms.eval_g_e     (0.5, 0.6, 0.7, 0.8, 3));
    OUTPUT(ms.eval_g_p     (0.5, 0.6, 0.7, 0.8, 1));
    OUTPUT(ms.eval_g_p     (0.5, 0.6, 0.7, 0.8, 2));
    OUTPUT(ms.eval_g_p     (0.5, 0.6, 0.7, 0.8, 3));
    OUTPUT(ms.eval_g_mu    (0.5, 0.6, 0.7, 0.8, 1));
    OUTPUT(ms.eval_g_mu    (0.5, 0.6, 0.7, 0.8, 2));
    OUTPUT(ms.eval_g_mu    (0.5, 0.6, 0.7, 0.8, 3));
    OUTPUT(ms.eval_q_rho   (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.eval_q_rho_u (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.eval_q_rho_v (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.eval_q_rho_w (0.5, 0.6, 0.7, 0.8));
    OUTPUT(ms.eval_q_rho_e (0.5, 0.6, 0.7, 0.8));

    return EXIT_SUCCESS;
}
