#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "nsctpl_fwd.hpp"
#include "nsctpl.hpp"
#include "test.hpp"

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

// Helper evaluate results against known good value.
// Use relative error (rather than absolute) as the values have large magnitudes
#define CHECK(expr,expected)                                                \
    std::cout << '\t' << #expr << " = " << (expr) << " has relerr "         \
              << ((expr)/(expected) - 1) << std::endl;

int main(int argc, char *argv[])
{
    std::cout.precision(std::numeric_limits<long double>::digits10);

    nsctpl::manufactured_solution<long double> ms;
    {
        prime_it<long double> p;
        ms.foreach_parameter(p);
    }

    std::cout << "Parameters after initialization:" << std::endl;
    ms.foreach_parameter(print_it<long double>);

    using nsctpl::test::x;
    using nsctpl::test::y;
    using nsctpl::test::z;
    using nsctpl::test::t;
    std::cout << "Evaluation takes place at" << std::endl;
    std::cout << "\tx = " << x << std::endl;
    std::cout << "\ty = " << y << std::endl;
    std::cout << "\tz = " << z << std::endl;
    std::cout << "\tt = " << t << std::endl;

    std::cout << "Analytic solution evaluations:" << std::endl;
    CHECK(ms.soln_rho    (x, y, z, t), nsctpl::test::rho);
    CHECK(ms.soln_rho._t (x, y, z, t), nsctpl::test::rho_t );
    CHECK(ms.soln_rho._x (x, y, z, t), nsctpl::test::rho_x );
    CHECK(ms.soln_rho._xx(x, y, z, t), nsctpl::test::rho_xx);
    CHECK(ms.soln_rho._xy(x, y, z, t), nsctpl::test::rho_xy);
    CHECK(ms.soln_rho._xz(x, y, z, t), nsctpl::test::rho_xz);
    CHECK(ms.soln_rho._y (x, y, z, t), nsctpl::test::rho_y );
    CHECK(ms.soln_rho._yy(x, y, z, t), nsctpl::test::rho_yy);
    CHECK(ms.soln_rho._yz(x, y, z, t), nsctpl::test::rho_yz);
    CHECK(ms.soln_rho._z (x, y, z, t), nsctpl::test::rho_z );
    CHECK(ms.soln_rho._zz(x, y, z, t), nsctpl::test::rho_zz);
    CHECK(ms.soln_u    (x, y, z, t), nsctpl::test::u   );
    CHECK(ms.soln_u._t (x, y, z, t), nsctpl::test::u_t );
    CHECK(ms.soln_u._x (x, y, z, t), nsctpl::test::u_x );
    CHECK(ms.soln_u._xx(x, y, z, t), nsctpl::test::u_xx);
    CHECK(ms.soln_u._xy(x, y, z, t), nsctpl::test::u_xy);
    CHECK(ms.soln_u._xz(x, y, z, t), nsctpl::test::u_xz);
    CHECK(ms.soln_u._y (x, y, z, t), nsctpl::test::u_y );
    CHECK(ms.soln_u._yy(x, y, z, t), nsctpl::test::u_yy);
    CHECK(ms.soln_u._yz(x, y, z, t), nsctpl::test::u_yz);
    CHECK(ms.soln_u._z (x, y, z, t), nsctpl::test::u_z );
    CHECK(ms.soln_u._zz(x, y, z, t), nsctpl::test::u_zz);
    CHECK(ms.soln_v    (x, y, z, t), nsctpl::test::v   );
    CHECK(ms.soln_v._t (x, y, z, t), nsctpl::test::v_t );
    CHECK(ms.soln_v._x (x, y, z, t), nsctpl::test::v_x );
    CHECK(ms.soln_v._xx(x, y, z, t), nsctpl::test::v_xx);
    CHECK(ms.soln_v._xy(x, y, z, t), nsctpl::test::v_xy);
    CHECK(ms.soln_v._xz(x, y, z, t), nsctpl::test::v_xz);
    CHECK(ms.soln_v._y (x, y, z, t), nsctpl::test::v_y );
    CHECK(ms.soln_v._yy(x, y, z, t), nsctpl::test::v_yy);
    CHECK(ms.soln_v._yz(x, y, z, t), nsctpl::test::v_yz);
    CHECK(ms.soln_v._z (x, y, z, t), nsctpl::test::v_z );
    CHECK(ms.soln_v._zz(x, y, z, t), nsctpl::test::v_zz);
    CHECK(ms.soln_w    (x, y, z, t), nsctpl::test::w   );
    CHECK(ms.soln_w._t (x, y, z, t), nsctpl::test::w_t );
    CHECK(ms.soln_w._x (x, y, z, t), nsctpl::test::w_x );
    CHECK(ms.soln_w._xx(x, y, z, t), nsctpl::test::w_xx);
    CHECK(ms.soln_w._xy(x, y, z, t), nsctpl::test::w_xy);
    CHECK(ms.soln_w._xz(x, y, z, t), nsctpl::test::w_xz);
    CHECK(ms.soln_w._y (x, y, z, t), nsctpl::test::w_y );
    CHECK(ms.soln_w._yy(x, y, z, t), nsctpl::test::w_yy);
    CHECK(ms.soln_w._yz(x, y, z, t), nsctpl::test::w_yz);
    CHECK(ms.soln_w._z (x, y, z, t), nsctpl::test::w_z );
    CHECK(ms.soln_w._zz(x, y, z, t), nsctpl::test::w_zz);
    CHECK(ms.soln_T    (x, y, z, t), nsctpl::test::T   );
    CHECK(ms.soln_T._t (x, y, z, t), nsctpl::test::T_t );
    CHECK(ms.soln_T._x (x, y, z, t), nsctpl::test::T_x );
    CHECK(ms.soln_T._xx(x, y, z, t), nsctpl::test::T_xx);
    CHECK(ms.soln_T._xy(x, y, z, t), nsctpl::test::T_xy);
    CHECK(ms.soln_T._xz(x, y, z, t), nsctpl::test::T_xz);
    CHECK(ms.soln_T._y (x, y, z, t), nsctpl::test::T_y );
    CHECK(ms.soln_T._yy(x, y, z, t), nsctpl::test::T_yy);
    CHECK(ms.soln_T._yz(x, y, z, t), nsctpl::test::T_yz);
    CHECK(ms.soln_T._z (x, y, z, t), nsctpl::test::T_z );
    CHECK(ms.soln_T._zz(x, y, z, t), nsctpl::test::T_zz);

    std::cout << "Analytically determined quantities" << std::endl;
    std::cout.precision(std::numeric_limits<long double>::digits10);
    CHECK(ms.eval_exact_rho(x, y, z, t), nsctpl::test::rho);
    CHECK(ms.eval_exact_u  (x, y, z, t), nsctpl::test::u);
    CHECK(ms.eval_exact_v  (x, y, z, t), nsctpl::test::v);
    CHECK(ms.eval_exact_w  (x, y, z, t), nsctpl::test::w);
    CHECK(ms.eval_exact_T  (x, y, z, t), nsctpl::test::T);
    CHECK(ms.eval_g_rho(x, y, z, t, 1),  nsctpl::test::rho_x);
    CHECK(ms.eval_g_rho(x, y, z, t, 2),  nsctpl::test::rho_y);
    CHECK(ms.eval_g_rho(x, y, z, t, 3),  nsctpl::test::rho_z);
    CHECK(ms.eval_g_u  (x, y, z, t, 1),  nsctpl::test::u_x);
    CHECK(ms.eval_g_u  (x, y, z, t, 2),  nsctpl::test::u_y);
    CHECK(ms.eval_g_u  (x, y, z, t, 3),  nsctpl::test::u_z);
    CHECK(ms.eval_g_v  (x, y, z, t, 1),  nsctpl::test::v_x);
    CHECK(ms.eval_g_v  (x, y, z, t, 2),  nsctpl::test::v_y);
    CHECK(ms.eval_g_v  (x, y, z, t, 3),  nsctpl::test::v_z);
    CHECK(ms.eval_g_w  (x, y, z, t, 1),  nsctpl::test::w_x);
    CHECK(ms.eval_g_w  (x, y, z, t, 2),  nsctpl::test::w_y);
    CHECK(ms.eval_g_w  (x, y, z, t, 3),  nsctpl::test::w_z);
    CHECK(ms.eval_g_T  (x, y, z, t, 1),  nsctpl::test::T_x);
    CHECK(ms.eval_g_T  (x, y, z, t, 2),  nsctpl::test::T_y);
    CHECK(ms.eval_g_T  (x, y, z, t, 3),  nsctpl::test::T_z);

    std::cout << "Quantities built from the analytical solutions" << std::endl;
    CHECK(ms.eval_exact_e (x, y, z, t), nsctpl::test::e);
    CHECK(ms.eval_exact_p (x, y, z, t), nsctpl::test::p);
    CHECK(ms.eval_exact_mu(x, y, z, t), nsctpl::test::mu);
    CHECK(ms.eval_g_e     (x, y, z, t, 1), nsctpl::test::e_x);
    CHECK(ms.eval_g_e     (x, y, z, t, 2), nsctpl::test::e_y);
    CHECK(ms.eval_g_e     (x, y, z, t, 3), nsctpl::test::e_z);
    CHECK(ms.eval_g_p     (x, y, z, t, 1), nsctpl::test::p_x);
    CHECK(ms.eval_g_p     (x, y, z, t, 2), nsctpl::test::p_y);
    CHECK(ms.eval_g_p     (x, y, z, t, 3), nsctpl::test::p_z);
    CHECK(ms.eval_g_mu    (x, y, z, t, 1), nsctpl::test::mu_x);
    CHECK(ms.eval_g_mu    (x, y, z, t, 2), nsctpl::test::mu_y);
    CHECK(ms.eval_g_mu    (x, y, z, t, 3), nsctpl::test::mu_z);
    CHECK(ms.eval_q_rho   (x, y, z, t), nsctpl::test::Q_rho);
    CHECK(ms.eval_q_rho_u (x, y, z, t), nsctpl::test::Q_rhou);
    CHECK(ms.eval_q_rho_v (x, y, z, t), nsctpl::test::Q_rhov);
    CHECK(ms.eval_q_rho_w (x, y, z, t), nsctpl::test::Q_rhow);
    CHECK(ms.eval_q_rho_e (x, y, z, t), nsctpl::test::Q_rhoe);

    return EXIT_SUCCESS;
}
