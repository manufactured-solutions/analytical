/**
 * POD-type providing operations for evaluating an analytical solution at a
 * given location and time.  The solution is of the form
 * \verbatim
 *       a_0                                         *sin(f_0 *t + g_0 )
 *     + a_x  * sin(b_x *x + c_x )                   *sin(f_x *t + g_x )
 *     + a_xy * sin(b_xy*x + c_xy)*sin(d_xy*y + e_xy)*sin(f_xy*t + g_xy)
 *     + a_xz * sin(b_xz*x + c_xz)*sin(d_xz*z + e_xz)*sin(f_xz*t + g_xz)
 *     + a_y  * sin(b_y *y + c_y )                   *sin(f_y *t + g_y )
 *     + a_yz * sin(b_yz*y + c_yz)*sin(d_yz*z + e_yz)*sin(f_yz*t + g_yz)
 *     + a_z  * sin(b_z *z + c_z )                   *sin(f_z *t + g_z )
 * \endverbatim
 *
 * \internal
 * Methods method names are non-traditional but align well with math notation.
 * Derivatives found using SymPy.  Linear storage within an internal array is
 * provided to simplify looping over all parameters within the solution.  Uses
 * X macros to simplify some repetitive coding
 * (http://drdobbs.com/blogs/cpp/228700289).
 */
template<typename FPT>
class soln {

public:

#define FOR_ALL_PARAMETERS(apply)                                                   \
    apply(a_,0),                                                                    \
    apply(a_,x), apply(a_,xy), apply(a_,xz), apply(a_,y), apply(a_,yz), apply(a_,z) \
    apply(b_,x), apply(b_,xy), apply(b_,xz), apply(b_,y), apply(b_,yz), apply(b_,z) \
    apply(c_,x), apply(c_,xy), apply(c_,xz), apply(c_,y), apply(c_,yz), apply(c_,z) \
                 apply(d_,xy), apply(d_,xz),              apply(d_,yz)              \
                 apply(e_,xy), apply(e_,xz),              apply(e_,yz)              \
    apply(f_,0),                                                                    \
    apply(f_,x), apply(f_,xy), apply(f_,xz), apply(f_,y), apply(f_,yz), apply(f_,z) \
    apply(g_,0),                                                                    \
    apply(g_,x), apply(g_,xy), apply(g_,xz), apply(g_,y), apply(g_,yz), apply(g_,z)

// Declare parameters
#define apply(pre,post) FPT pre##post;
    FOR_ALL_PARAMETERS(apply)
#undef apply


    FPT a_0, a_x, a_xy, a_xz, a_y, a_yz, a_z;
    FPT      b_x, b_xy, b_xz, b_y, b_yz, b_z;
    FPT      c_x, c_xy, c_xz, c_y, c_yz, c_z;
    FPT           d_xy, d_xz,      d_yz     ;
    FPT           e_xy, e_xz,      e_yz     ;
    FPT f_0, f_x, f_xy, f_xz, f_y, f_yz, f_z;
    FPT g_0, g_x, g_xy, g_xz, g_y, g_yz, g_z;

    //! Constructor which NaNs parameters to prevent non-initialization bugs
    soln() :
#define apply(pre,post) pre##post(std::numeric_limits<FPT>::quiet_NaN()),
    FOR_ALL_PARAMETERS(apply)
#undef apply
    {}

    //! Evaluate the solution
    FPT operator()(const FPT x, const FPT y, const FPT z, const FPT t);

    //! Evaluate the solution's derivative with respect to time
    FPT _t(const FPT x, const FPT y, const FPT z, const FPT t);

    //! Evaluate the solution's derivative with respect to \c x
    FPT _x(const FPT x, const FPT y, const FPT z, const FPT t);

    //! Evaluate the solution's second derivative with respect to \c x
    FPT _xx(const FPT x, const FPT y, const FPT z, const FPT t);

    //! Evaluate the solution's derivative with respect to \c x and \c y
    FPT _xy(const FPT x, const FPT y, const FPT z, const FPT t);

    //! Evaluate the solution's derivative with respect to \c x and \c z
    FPT _xz(const FPT x, const FPT y, const FPT z, const FPT t);

    //! Evaluate the solution's derivative with respect to \c y
    FPT _y(const FPT x, const FPT y, const FPT z, const FPT t);

    //! Evaluate the solution's second derivative with respect to \c y
    FPT _yy(const FPT x, const FPT y, const FPT z, const FPT t);

    //! Evaluate the solution's derivative with respect to \c y and \c z
    FPT _yz(const FPT x, const FPT y, const FPT z, const FPT t);

    //! Evaluate the solution's derivative with respect to \c z
    FPT _z(const FPT x, const FPT y, const FPT z, const FPT t);

    //! Evaluate the solution's second derivative with respect to \c z
    FPT _zz(const FPT x, const FPT y, const FPT z, const FPT t);

    //! Evaluate solution and all available derivatives
    void evaluate(const FPT x, const FPT y, const FPT z, const FPT t,
                  FPT& phi, FPT& phi_t,
                  FPT& phi_x, FPT& phi_xx, FPT& phi_xy, FPT& phi_xz,
                  FPT& phi_y, FPT& phi_yy, FPT& phi_yz,
                  FPT& phi_z, FPT& phi_zz) {

        phi    = this->operator()(x, y, z, t);
        phi_t  = this->_t        (x, y, z, t);
        phi_x  = this->_x        (x, y, z, t);
        phi_xx = this->_xx       (x, y, z, t);
        phi_xy = this->_xy       (x, y, z, t);
        phi_xz = this->_xz       (x, y, z, t);
        phi_y  = this->_y        (x, y, z, t);
        phi_yy = this->_yy       (x, y, z, t);
        phi_yz = this->_yz       (x, y, z, t);
        phi_z  = this->_z        (x, y, z, t);
        phi_zz = this->_zz       (x, y, z, t);
    }

// Avoid polluting the global namespace
#undef FOR_ALL_PARAMETERS

}; // end class

template<typename FPT>
FPT soln<FPT>::operator()(const FPT x, const FPT y, const FPT z, const FPT t)
{
    return a_0*std::sin(g_0 + f_0*t) + a_x*std::sin(c_x + b_x*x)*std::sin(g_x + f_x*t) + a_y*std::sin(g_y + f_y*t)*std::sin(c_y + b_y*y) + a_z*std::sin(g_z + f_z*t)*std::sin(c_z + b_z*z) + a_xy*std::sin(g_xy + f_xy*t)*std::sin(c_xy + b_xy*x)*std::sin(e_xy + d_xy*y) + a_xz*std::sin(c_xz + b_xz*x)*std::sin(g_xz + f_xz*t)*std::sin(e_xz + d_xz*z) + a_yz*std::sin(c_yz + b_yz*y)*std::sin(e_yz + d_yz*z)*std::sin(g_yz + f_yz*t);
}

template<typename FPT>
FPT soln<FPT>::_t(const FPT x, const FPT y, const FPT z, const FPT t)
{
    return a_0*f_0*std::cos(g_0 + f_0*t) + a_x*f_x*std::cos(g_x + f_x*t)*std::sin(c_x + b_x*x) + a_y*f_y*std::cos(g_y + f_y*t)*std::sin(c_y + b_y*y) + a_z*f_z*std::cos(g_z + f_z*t)*std::sin(c_z + b_z*z) + a_xy*f_xy*std::cos(g_xy + f_xy*t)*std::sin(c_xy + b_xy*x)*std::sin(e_xy + d_xy*y) + a_xz*f_xz*std::cos(g_xz + f_xz*t)*std::sin(c_xz + b_xz*x)*std::sin(e_xz + d_xz*z) + a_yz*f_yz*std::cos(g_yz + f_yz*t)*std::sin(c_yz + b_yz*y)*std::sin(e_yz + d_yz*z);
}

template<typename FPT>
FPT soln<FPT>::_x(const FPT x, const FPT y, const FPT z, const FPT t)
{
    return a_x*b_x*std::cos(c_x + b_x*x)*std::sin(g_x + f_x*t) + a_xy*b_xy*std::cos(c_xy + b_xy*x)*std::sin(g_xy + f_xy*t)*std::sin(e_xy + d_xy*y) + a_xz*b_xz*std::cos(c_xz + b_xz*x)*std::sin(g_xz + f_xz*t)*std::sin(e_xz + d_xz*z);
}

template<typename FPT>
FPT soln<FPT>::_xx(const FPT x, const FPT y, const FPT z, const FPT t)
{
    return -a_x*std::pow(b_x,2)*std::sin(c_x + b_x*x)*std::sin(g_x + f_x*t) - a_xy*std::pow(b_xy,2)*std::sin(g_xy + f_xy*t)*std::sin(c_xy + b_xy*x)*std::sin(e_xy + d_xy*y) - a_xz*std::pow(b_xz,2)*std::sin(c_xz + b_xz*x)*std::sin(g_xz + f_xz*t)*std::sin(e_xz + d_xz*z);
}

template<typename FPT>
FPT soln<FPT>::_xy(const FPT x, const FPT y, const FPT z, const FPT t)
{
    return a_xy*b_xy*d_xy*std::cos(c_xy + b_xy*x)*std::cos(e_xy + d_xy*y)*std::sin(g_xy + f_xy*t);
}

template<typename FPT>
FPT soln<FPT>::_xz(const FPT x, const FPT y, const FPT z, const FPT t)
{
    return a_xz*b_xz*d_xz*std::cos(c_xz + b_xz*x)*std::cos(e_xz + d_xz*z)*std::sin(g_xz + f_xz*t);
}

template<typename FPT>
FPT soln<FPT>::_y(const FPT x, const FPT y, const FPT z, const FPT t)
{
    return a_y*b_y*std::cos(c_y + b_y*y)*std::sin(g_y + f_y*t) + a_xy*d_xy*std::cos(e_xy + d_xy*y)*std::sin(g_xy + f_xy*t)*std::sin(c_xy + b_xy*x) + a_yz*b_yz*std::cos(c_yz + b_yz*y)*std::sin(e_yz + d_yz*z)*std::sin(g_yz + f_yz*t);
}

template<typename FPT>
FPT soln<FPT>::_yy(const FPT x, const FPT y, const FPT z, const FPT t)
{
    return -a_y*std::pow(b_y,2)*std::sin(g_y + f_y*t)*std::sin(c_y + b_y*y) - a_xy*std::pow(d_xy,2)*std::sin(g_xy + f_xy*t)*std::sin(c_xy + b_xy*x)*std::sin(e_xy + d_xy*y) - a_yz*std::pow(b_yz,2)*std::sin(c_yz + b_yz*y)*std::sin(e_yz + d_yz*z)*std::sin(g_yz + f_yz*t);
}

template<typename FPT>
FPT soln<FPT>::_yz(const FPT x, const FPT y, const FPT z, const FPT t)
{
    return a_yz*b_yz*d_yz*std::cos(c_yz + b_yz*y)*std::cos(e_yz + d_yz*z)*std::sin(g_yz + f_yz*t);
}

template<typename FPT>
FPT soln<FPT>::_z(const FPT x, const FPT y, const FPT z, const FPT t)
{
    return a_z*b_z*std::cos(c_z + b_z*z)*std::sin(g_z + f_z*t) + a_xz*d_xz*std::cos(e_xz + d_xz*z)*std::sin(c_xz + b_xz*x)*std::sin(g_xz + f_xz*t) + a_yz*d_yz*std::cos(e_yz + d_yz*z)*std::sin(c_yz + b_yz*y)*std::sin(g_yz + f_yz*t);
}

template<typename FPT>
FPT soln<FPT>::_zz(const FPT x, const FPT y, const FPT z, const FPT t)
{
    return -a_z*std::pow(b_z,2)*std::sin(g_z + f_z*t)*std::sin(c_z + b_z*z) - a_xz*std::pow(d_xz,2)*std::sin(c_xz + b_xz*x)*std::sin(g_xz + f_xz*t)*std::sin(e_xz + d_xz*z) - a_yz*std::pow(d_yz,2)*std::sin(c_yz + b_yz*y)*std::sin(e_yz + d_yz*z)*std::sin(g_yz + f_yz*t);
}
