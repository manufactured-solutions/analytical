 /* SWIG declarations for building nsctpl module     */
 /* See SWIG documentation at http://www.swig.org/   */
 /* Especially the SWIG and C++ details on templates */

%module nsctpl

%{
#define SWIG_FILE_WITH_INIT
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include "nsctpl.hpp"
%}

%include "nsctpl_fwd.hpp"

namespace nsctpl {

%extend primitive {

    %template(__call__) operator()<double,double,double,double>;
    %template(_t      ) _t        <double,double,double,double>;
    %template(_x      ) _x        <double,double,double,double>;
    %template(_xx     ) _xx       <double,double,double,double>;
    %template(_xy     ) _xy       <double,double,double,double>;
    %template(_xz     ) _xz       <double,double,double,double>;
    %template(_y      ) _y        <double,double,double,double>;
    %template(_yy     ) _yy       <double,double,double,double>;
    %template(_yz     ) _yz       <double,double,double,double>;
    %template(_z      ) _z        <double,double,double,double>;
    %template(_zz     ) _zz       <double,double,double,double>;
}

%template(primitive_double) primitive<double>;
%pythoncode %{
    primitive = primitive_double
%}

%extend manufactured_solution {
    %template(grad_rho) grad_rho<double,double,double,double>;
    %template(grad_u  ) grad_u  <double,double,double,double>;
    %template(grad_v  ) grad_v  <double,double,double,double>;
    %template(grad_w  ) grad_w  <double,double,double,double>;
    %template(grad_T  ) grad_T  <double,double,double,double>;
    %template(e       ) e       <double,double,double,double>;
    %template(p       ) p       <double,double,double,double>;
    %template(mu      ) mu      <double,double,double,double>;
    %template(rhou    ) rhou    <double,double,double,double>;
    %template(rhov    ) rhov    <double,double,double,double>;
    %template(rhow    ) rhow    <double,double,double,double>;
    %template(rhoe    ) rhoe    <double,double,double,double>;
    %template(grad_e  ) grad_e  <double,double,double,double>;
    %template(grad_p  ) grad_p  <double,double,double,double>;
    %template(grad_mu ) grad_mu <double,double,double,double>;
    %template(Q_rho   ) Q_rho   <double,double,double,double>;
    %template(Q_rhou  ) Q_rhou  <double,double,double,double>;
    %template(Q_rhov  ) Q_rhov  <double,double,double,double>;
    %template(Q_rhow  ) Q_rhow  <double,double,double,double>;
    %template(Q_rhoe  ) Q_rhoe  <double,double,double,double>;
}

%template(manufactured_solution_double) manufactured_solution<double>;
%pythoncode %{
    manufactured_solution = manufactured_solution_double
%}

}
