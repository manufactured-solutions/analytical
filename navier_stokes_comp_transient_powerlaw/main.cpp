#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "nsctpl.hpp"

// Explicitly instantiate the solution for several floating point types
template class solution<float>;
template class solution<double>;
template class solution<long double>;

// Used to output each solution parameter on std::cout
void printer(const std::string &name, double value) {
    std::cout << '\t' << name << " = " << value << std::endl;
}

// Set each parameter value to be a (very contribed) hash of its name
void initter(const std::string &name, double& value) {
    value = 0;
    for (std::size_t i = 0; i < name.length(); ++i) {
        value += name[i];
    }
    value /= 1000;
}

int main(int argc, char *argv[])
{
    const double x = 0.5, y = 0.6, z = 0.7, t = 0.8;
    solution<double> soln_rho("rho");

    std::cout << "Parameters after construction:" << std::endl;
    soln_rho.forall_parameters(printer);

    std::cout << "Evaluation after construction:" << std::endl;
    std::cout << "\trho    = " << soln_rho    (x, y, z, t) << std::endl;
    std::cout << "\trho_t  = " << soln_rho._t (x, y, z, t) << std::endl;
    std::cout << "\trho_x  = " << soln_rho._x (x, y, z, t) << std::endl;
    std::cout << "\trho_xx = " << soln_rho._xx(x, y, z, t) << std::endl;
    std::cout << "\trho_xy = " << soln_rho._xy(x, y, z, t) << std::endl;
    std::cout << "\trho_xz = " << soln_rho._xz(x, y, z, t) << std::endl;
    std::cout << "\trho_y  = " << soln_rho._y (x, y, z, t) << std::endl;
    std::cout << "\trho_yy = " << soln_rho._yy(x, y, z, t) << std::endl;
    std::cout << "\trho_yz = " << soln_rho._yz(x, y, z, t) << std::endl;
    std::cout << "\trho_z  = " << soln_rho._z (x, y, z, t) << std::endl;
    std::cout << "\trho_zz = " << soln_rho._zz(x, y, z, t) << std::endl;

    std::cout << "Parameters after initializing:" << std::endl;
    soln_rho.forall_parameters(initter);
    soln_rho.a_0  = 11;
    soln_rho.b_x  = 13;
    soln_rho.c_xy = 17;
    soln_rho.d_xz = 19;
    soln_rho.e_yz = 23;
    soln_rho.f_z  = 29;
    soln_rho.g_z  = 31;
    soln_rho.forall_parameters(printer);

    std::cout << "Evaluation after initializing:" << std::endl;
    std::cout << "\trho    = " << soln_rho    (x, y, z, t) << std::endl;
    std::cout << "\trho_t  = " << soln_rho._t (x, y, z, t) << std::endl;
    std::cout << "\trho_x  = " << soln_rho._x (x, y, z, t) << std::endl;
    std::cout << "\trho_xx = " << soln_rho._xx(x, y, z, t) << std::endl;
    std::cout << "\trho_xy = " << soln_rho._xy(x, y, z, t) << std::endl;
    std::cout << "\trho_xz = " << soln_rho._xz(x, y, z, t) << std::endl;
    std::cout << "\trho_y  = " << soln_rho._y (x, y, z, t) << std::endl;
    std::cout << "\trho_yy = " << soln_rho._yy(x, y, z, t) << std::endl;
    std::cout << "\trho_yz = " << soln_rho._yz(x, y, z, t) << std::endl;
    std::cout << "\trho_z  = " << soln_rho._z (x, y, z, t) << std::endl;
    std::cout << "\trho_zz = " << soln_rho._zz(x, y, z, t) << std::endl;

    return EXIT_SUCCESS;
}
