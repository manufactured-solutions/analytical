#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <typeinfo>
#include <tr1/memory>

#include "nsctpl_rholut.hpp"
#include "test.hpp"

#define STRINGIFY(foo) #foo

// Explicitly instantiate the primitive solution to ensure compilation
template class nsctpl_rholut::primitive<float>;
template class nsctpl_rholut::primitive<double>;
template class nsctpl_rholut::primitive<long double>;

// Explicitly instantiate the manufactured forcing to ensure compilation
template class nsctpl_rholut::manufactured_solution<float,1>;
template class nsctpl_rholut::manufactured_solution<double,1>;
template class nsctpl_rholut::manufactured_solution<long double,1>;

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

    prime_it() : i(new std::size_t) { *i = 0; }

    void operator()(const std::string&, Scalar& value) {
        assert(*i < sizeof(primes)/sizeof(primes[0]));
        value = primes[(*i)++];
    }

private:
    std::tr1::shared_ptr<std::size_t> i; // shared_ptr so copies share state
};

// Helper evaluate results against known good value.
// Use relative error (rather than absolute) as the values have large magnitudes
#define CHECK(expr,expected)                                                \
    std::cout << '\t' << #expr << " = " << (expr) << " has relerr "         \
              << ((expr)/(expected) - 1) << std::endl;

template<typename Scalar>
void run_tests()
{
    std::cout.precision(std::numeric_limits<Scalar>::digits10);

    nsctpl_rholut::manufactured_solution<Scalar,1> ms;
    {
        prime_it<Scalar> p;
        ms.foreach_parameter(p);
    }

    using nsctpl_rholut::test::x;
    using nsctpl_rholut::test::y;
    using nsctpl_rholut::test::z;
    using nsctpl_rholut::test::t;

    std::cout << std::endl
              << "Analytic solution evaluations for type "
              << typeid(Scalar).name()
              << " (digits10 = " << std::numeric_limits<Scalar>::digits10 << ")"
              << std::endl;
    CHECK(ms.rho    (x, y, z, t), nsctpl_rholut::test::rho);
    CHECK(ms.rho._t (x, y, z, t), nsctpl_rholut::test::rho_t );
    CHECK(ms.rho._x (x, y, z, t), nsctpl_rholut::test::rho_x );
    CHECK(ms.rho._xx(x, y, z, t), nsctpl_rholut::test::rho_xx);
    CHECK(ms.rho._xy(x, y, z, t), nsctpl_rholut::test::rho_xy);
    CHECK(ms.rho._xz(x, y, z, t), nsctpl_rholut::test::rho_xz);
    CHECK(ms.rho._y (x, y, z, t), nsctpl_rholut::test::rho_y );
    CHECK(ms.rho._yy(x, y, z, t), nsctpl_rholut::test::rho_yy);
    CHECK(ms.rho._yz(x, y, z, t), nsctpl_rholut::test::rho_yz);
    CHECK(ms.rho._z (x, y, z, t), nsctpl_rholut::test::rho_z );
    CHECK(ms.rho._zz(x, y, z, t), nsctpl_rholut::test::rho_zz);
    CHECK(ms.u    (x, y, z, t), nsctpl_rholut::test::u   );
    CHECK(ms.u._t (x, y, z, t), nsctpl_rholut::test::u_t );
    CHECK(ms.u._x (x, y, z, t), nsctpl_rholut::test::u_x );
    CHECK(ms.u._xx(x, y, z, t), nsctpl_rholut::test::u_xx);
    CHECK(ms.u._xy(x, y, z, t), nsctpl_rholut::test::u_xy);
    CHECK(ms.u._xz(x, y, z, t), nsctpl_rholut::test::u_xz);
    CHECK(ms.u._y (x, y, z, t), nsctpl_rholut::test::u_y );
    CHECK(ms.u._yy(x, y, z, t), nsctpl_rholut::test::u_yy);
    CHECK(ms.u._yz(x, y, z, t), nsctpl_rholut::test::u_yz);
    CHECK(ms.u._z (x, y, z, t), nsctpl_rholut::test::u_z );
    CHECK(ms.u._zz(x, y, z, t), nsctpl_rholut::test::u_zz);
    CHECK(ms.v    (x, y, z, t), nsctpl_rholut::test::v   );
    CHECK(ms.v._t (x, y, z, t), nsctpl_rholut::test::v_t );
    CHECK(ms.v._x (x, y, z, t), nsctpl_rholut::test::v_x );
    CHECK(ms.v._xx(x, y, z, t), nsctpl_rholut::test::v_xx);
    CHECK(ms.v._xy(x, y, z, t), nsctpl_rholut::test::v_xy);
    CHECK(ms.v._xz(x, y, z, t), nsctpl_rholut::test::v_xz);
    CHECK(ms.v._y (x, y, z, t), nsctpl_rholut::test::v_y );
    CHECK(ms.v._yy(x, y, z, t), nsctpl_rholut::test::v_yy);
    CHECK(ms.v._yz(x, y, z, t), nsctpl_rholut::test::v_yz);
    CHECK(ms.v._z (x, y, z, t), nsctpl_rholut::test::v_z );
    CHECK(ms.v._zz(x, y, z, t), nsctpl_rholut::test::v_zz);
    CHECK(ms.w    (x, y, z, t), nsctpl_rholut::test::w   );
    CHECK(ms.w._t (x, y, z, t), nsctpl_rholut::test::w_t );
    CHECK(ms.w._x (x, y, z, t), nsctpl_rholut::test::w_x );
    CHECK(ms.w._xx(x, y, z, t), nsctpl_rholut::test::w_xx);
    CHECK(ms.w._xy(x, y, z, t), nsctpl_rholut::test::w_xy);
    CHECK(ms.w._xz(x, y, z, t), nsctpl_rholut::test::w_xz);
    CHECK(ms.w._y (x, y, z, t), nsctpl_rholut::test::w_y );
    CHECK(ms.w._yy(x, y, z, t), nsctpl_rholut::test::w_yy);
    CHECK(ms.w._yz(x, y, z, t), nsctpl_rholut::test::w_yz);
    CHECK(ms.w._z (x, y, z, t), nsctpl_rholut::test::w_z );
    CHECK(ms.w._zz(x, y, z, t), nsctpl_rholut::test::w_zz);
    CHECK(ms.T    (x, y, z, t), nsctpl_rholut::test::T   );
    CHECK(ms.T._t (x, y, z, t), nsctpl_rholut::test::T_t );
    CHECK(ms.T._x (x, y, z, t), nsctpl_rholut::test::T_x );
    CHECK(ms.T._xx(x, y, z, t), nsctpl_rholut::test::T_xx);
    CHECK(ms.T._xy(x, y, z, t), nsctpl_rholut::test::T_xy);
    CHECK(ms.T._xz(x, y, z, t), nsctpl_rholut::test::T_xz);
    CHECK(ms.T._y (x, y, z, t), nsctpl_rholut::test::T_y );
    CHECK(ms.T._yy(x, y, z, t), nsctpl_rholut::test::T_yy);
    CHECK(ms.T._yz(x, y, z, t), nsctpl_rholut::test::T_yz);
    CHECK(ms.T._z (x, y, z, t), nsctpl_rholut::test::T_z );
    CHECK(ms.T._zz(x, y, z, t), nsctpl_rholut::test::T_zz);

    std::cout << std::endl
              << "Analytically determined quantities for "
              << typeid(Scalar).name()
              << " (digits10 = " << std::numeric_limits<Scalar>::digits10 << ")"
              << std::endl;
    CHECK(ms.grad_rho(x, y, z, t, 1), nsctpl_rholut::test::rho_x);
    CHECK(ms.grad_rho(x, y, z, t, 2), nsctpl_rholut::test::rho_y);
    CHECK(ms.grad_rho(x, y, z, t, 3), nsctpl_rholut::test::rho_z);
    CHECK(ms.grad_u  (x, y, z, t, 1), nsctpl_rholut::test::u_x);
    CHECK(ms.grad_u  (x, y, z, t, 2), nsctpl_rholut::test::u_y);
    CHECK(ms.grad_u  (x, y, z, t, 3), nsctpl_rholut::test::u_z);
    CHECK(ms.grad_v  (x, y, z, t, 1), nsctpl_rholut::test::v_x);
    CHECK(ms.grad_v  (x, y, z, t, 2), nsctpl_rholut::test::v_y);
    CHECK(ms.grad_v  (x, y, z, t, 3), nsctpl_rholut::test::v_z);
    CHECK(ms.grad_w  (x, y, z, t, 1), nsctpl_rholut::test::w_x);
    CHECK(ms.grad_w  (x, y, z, t, 2), nsctpl_rholut::test::w_y);
    CHECK(ms.grad_w  (x, y, z, t, 3), nsctpl_rholut::test::w_z);
    CHECK(ms.grad_T  (x, y, z, t, 1), nsctpl_rholut::test::T_x);
    CHECK(ms.grad_T  (x, y, z, t, 2), nsctpl_rholut::test::T_y);
    CHECK(ms.grad_T  (x, y, z, t, 3), nsctpl_rholut::test::T_z);

    std::cout << std::endl
              << "Quantities built from the analytical solutions for "
              << typeid(Scalar).name()
              << " (digits10 = " << std::numeric_limits<Scalar>::digits10 << ")"
              << std::endl;
    CHECK(ms.e      (x, y, z, t),    nsctpl_rholut::test::e);
    CHECK(ms.p      (x, y, z, t),    nsctpl_rholut::test::p);
    CHECK(ms.mu     (x, y, z, t),    nsctpl_rholut::test::mu);
    CHECK(ms.rhou   (x, y, z, t),    nsctpl_rholut::test::rhou);
    CHECK(ms.rhov   (x, y, z, t),    nsctpl_rholut::test::rhov);
    CHECK(ms.rhow   (x, y, z, t),    nsctpl_rholut::test::rhow);
    CHECK(ms.rhoe   (x, y, z, t),    nsctpl_rholut::test::rhoe);
    CHECK(ms.grad_e (x, y, z, t, 1), nsctpl_rholut::test::e_x);
    CHECK(ms.grad_e (x, y, z, t, 2), nsctpl_rholut::test::e_y);
    CHECK(ms.grad_e (x, y, z, t, 3), nsctpl_rholut::test::e_z);
    CHECK(ms.grad_p (x, y, z, t, 1), nsctpl_rholut::test::p_x);
    CHECK(ms.grad_p (x, y, z, t, 2), nsctpl_rholut::test::p_y);
    CHECK(ms.grad_p (x, y, z, t, 3), nsctpl_rholut::test::p_z);
    CHECK(ms.grad_mu(x, y, z, t, 1), nsctpl_rholut::test::mu_x);
    CHECK(ms.grad_mu(x, y, z, t, 2), nsctpl_rholut::test::mu_y);
    CHECK(ms.grad_mu(x, y, z, t, 3), nsctpl_rholut::test::mu_z);
    CHECK(ms.Q_rho  (x, y, z, t),    nsctpl_rholut::test::Q_rho);
    CHECK(ms.Q_rhou (x, y, z, t),    nsctpl_rholut::test::Q_rhou);
    CHECK(ms.Q_rhov (x, y, z, t),    nsctpl_rholut::test::Q_rhov);
    CHECK(ms.Q_rhow (x, y, z, t),    nsctpl_rholut::test::Q_rhow);
    CHECK(ms.Q_rhoe (x, y, z, t),    nsctpl_rholut::test::Q_rhoe);

    // Test Q_conservative results: simultaneous Q_rho{,u,v,w,e} computation
    Scalar Q_conservative_rho  = 0;
    Scalar Q_conservative_rhou = 0;
    Scalar Q_conservative_rhov = 0;
    Scalar Q_conservative_rhow = 0;
    Scalar Q_conservative_rhoe = 0;
    ms.Q_conservative(x, y, z, t, Q_conservative_rho,
                                  Q_conservative_rhou,
                                  Q_conservative_rhov,
                                  Q_conservative_rhow,
                                  Q_conservative_rhoe);
    CHECK(Q_conservative_rho , nsctpl_rholut::test::Q_rho);
    CHECK(Q_conservative_rhou, nsctpl_rholut::test::Q_rhou);
    CHECK(Q_conservative_rhov, nsctpl_rholut::test::Q_rhov);
    CHECK(Q_conservative_rhow, nsctpl_rholut::test::Q_rhow);
    CHECK(Q_conservative_rhoe, nsctpl_rholut::test::Q_rhoe);
}

int main(int argc, char *argv[])
{
    using std::cout;
    using std::endl;

    cout << "navier_stokes_comp_transient_powerlaw_rholut implementation test" << endl
         << '\t' << STRINGIFY($Id$) << endl;

    cout << endl
         << "Solution parameters after initialization:" << endl;
    nsctpl_rholut::manufactured_solution<double> ms;
    prime_it<double> p;
    ms.foreach_parameter(p);
    ms.foreach_parameter(print_it<double>);

    cout << endl
         << "Evaluation takes place at" << endl
         << "\tx = " << nsctpl_rholut::test::x << endl
         << "\ty = " << nsctpl_rholut::test::y << endl
         << "\tz = " << nsctpl_rholut::test::z << endl
         << "\tt = " << nsctpl_rholut::test::t << endl;

    run_tests<float>();
    run_tests<double>();
    run_tests<long double>();

    return EXIT_SUCCESS;
}
