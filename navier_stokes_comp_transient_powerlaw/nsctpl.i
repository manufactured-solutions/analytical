// SWIG declarations for building nsctpl module
// See SWIG documentation at http://www.swig.org/
// Especially the SWIG and C++ details on templates

%module nsctpl

// Need double* methods to interact with nsctpl::primitive Scalar* members
%include "cpointer.i"
%pointer_functions(double, doublep);

// Verbatim #includes passed into C++ compilation taken from nsctpl_fwd.hpp
%{
#define SWIG_FILE_WITH_INIT
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include "nsctpl.hpp"
%}

// Verbatim definitions for within the Python module declaration
%pythoncode %{

    # A metaclass that adds enthought.traits.api.HasTraits as an ancestor
    # See http://onlamp.com/lpt/a/3388 for an intro to Python metaclasses and
    # http://code.activestate.com/recipes/204197-solving-the-metaclass-conflict/
    from enthought.traits.api import HasTraits
    from enthought.traits.has_traits import MetaHasTraits
    class AddHasTraitsAncestor(MetaHasTraits):
        def __new__(cls, name, bases, dct):
            return MetaHasTraits.__new__(cls, name, cls.modify(bases), dct)

        def __init__(cls, name, bases, dct):
            super(AddHasTraitsAncestor, cls).__init__(name, cls.modify(bases), dct)

        @staticmethod
        def modify(bases):
            if (bases == (object,)):
                return (HasTraits,)
            else:
                new_bases = list(bases)
                new_bases.insert(0, HasTraits)
                return tuple(new_bases)

%}

// Have SWIG parse nsctpl namespace forward declarations
%include "nsctpl_fwd.hpp"

namespace nsctpl {

// Instantiate templated nsctpl::primitive members for doubles
%extend primitive {

    // Cause Python class to inherit from enthought.traits.api.HasTraits
    %pythoncode %{
        __metaclass__ = AddHasTraitsAncestor
    %}

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

// Expose template instantiation of nsctpl::primitive for doubles
%template(primitive_double) primitive<double>;
%pythoncode %{
    primitive = primitive_double
%}

// Instantiate templated nsctpl::manufactured_solution members for doubles
%extend manufactured_solution {

    // Cause Python class to inherit from enthought.traits.api.HasTraits
    %pythoncode %{
        __metaclass__ = AddHasTraitsAncestor
    %}

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

// Expose template instantiation of nsctpl::manufactured_solution for doubles
%template(manufactured_solution_double) manufactured_solution<double>;
%pythoncode %{
    manufactured_solution = manufactured_solution_double
%}

}
