/* SWIG declarations for building nsctpl module   */
/* See SWIG documentation at http://www.swig.org/ */

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
    %template()   primitive_solution<double>;
    %template()   generic_manufactured_solution<primitive_solution,double,1>;
    %template(ms) manufactured_solution<double>;
}
