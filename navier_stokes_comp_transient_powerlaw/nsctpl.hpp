// Using this file will require #include-ing <cmath>, <sstream>, and <string>.
// These headers are not included here to allow #include-ing this file inside
// any namespace, as is the consistent use of ::std instead of merely std.

/**
 * POD-type providing operations for evaluating an analytical solution at a
 * given location and time.  The solution is of the form
 * \verbatim
 *       a_0                                         *cos(f_0 *t + g_0 )
 *     + a_x  * cos(b_x *x + c_x )                   *cos(f_x *t + g_x )
 *     + a_xy * cos(b_xy*x + c_xy)*cos(d_xy*y + e_xy)*cos(f_xy*t + g_xy)
 *     + a_xz * cos(b_xz*x + c_xz)*cos(d_xz*z + e_xz)*cos(f_xz*t + g_xz)
 *     + a_y  * cos(b_y *y + c_y )                   *cos(f_y *t + g_y )
 *     + a_yz * cos(b_yz*y + c_yz)*cos(d_yz*z + e_yz)*cos(f_yz*t + g_yz)
 *     + a_z  * cos(b_z *z + c_z )                   *cos(f_z *t + g_z )
 * \endverbatim
 * The member method names are non-traditional but permitted by the language
 * standards and well-aligned with mathematical notation.  The
 * forall_parameters() member allows invoking an operation on each solution
 * parameter to aid parameter registration, initialization, or output.
 * Parameters may have an infix name added for use with forall_parameters().
 * For example, setting the name 'phi' will cause forall_parameters() to report
 * names like 'a_phix'.
 */
template<typename Scalar>
class solution {

public:

     // X macro per http://drdobbs.com/blogs/cpp/228700289 which will
     // apply a macro on all solution parameter prefix/suffix pairs.
#define FOR_ALL_SOLUTION_PARAMETERS(apply)                                     \
    apply(a_,0)                                                                \
    apply(a_,x) apply(a_,xy) apply(a_,xz) apply(a_,y) apply(a_,yz) apply(a_,z) \
    apply(b_,x) apply(b_,xy) apply(b_,xz) apply(b_,y) apply(b_,yz) apply(b_,z) \
    apply(c_,x) apply(c_,xy) apply(c_,xz) apply(c_,y) apply(c_,yz) apply(c_,z) \
                apply(d_,xy) apply(d_,xz)             apply(d_,yz)             \
                apply(e_,xy) apply(e_,xz)             apply(e_,yz)             \
    apply(f_,0)                                                                \
    apply(f_,x) apply(f_,xy) apply(f_,xz) apply(f_,y) apply(f_,yz) apply(f_,z) \
    apply(g_,0)                                                                \
    apply(g_,x) apply(g_,xy) apply(g_,xz) apply(g_,y) apply(g_,yz) apply(g_,z)

    // Declare all solution parameters as members, e.g. 'a_xy'
#define APPLY(pre,suf) Scalar pre##suf;
    FOR_ALL_SOLUTION_PARAMETERS(APPLY)
#undef APPLY

    //! The name used infix in forall_parameters' names.  For example,
    //! name = 'phi' implies parameters names like 'a_phix'.
    const ::std::string name;

    //! Construct an instance using \c name in the reported parameter names.
    //! All parameters set to zero at construction time.
    explicit solution(const ::std::string &name = "")
#define APPLY(pre,suf) pre##suf(static_cast<Scalar>(0)),
        : FOR_ALL_SOLUTION_PARAMETERS(APPLY)  // has trailing comma
#undef APPLY
          name(name)
    {}

#define APPLY_STRINGIFY(s) #s
#define APPLY(pre,suf)                                                    \
        os.clear(); os.str("");                                           \
        os << APPLY_STRINGIFY(pre) << this->name << APPLY_STRINGIFY(suf); \
        f(os.str(), this->pre##suf);

    //! Invoke the binary function f on each parameter name and its value.
    template<typename BinaryFunction>
    void forall_parameters(BinaryFunction f) const {
        ::std::ostringstream os;
        FOR_ALL_SOLUTION_PARAMETERS(APPLY)
    }

    //! Invoke the binary function f on each parameter name and its value.
    template<typename BinaryFunction>
    void forall_parameters(BinaryFunction f) {
        ::std::ostringstream os;
        FOR_ALL_SOLUTION_PARAMETERS(APPLY)
    }

#undef APPLY
#undef APPLY_STRINGIFY

// Avoid polluting the global namespace with an implementation detail
#undef FOR_ALL_SOLUTION_PARAMETERS

    //! Evaluate the solution
    Scalar operator()(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const;

    //! Evaluate the solution's derivative with respect to time
    Scalar _t(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const;

    //! Evaluate the solution's derivative with respect to \c x
    Scalar _x(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const;

    //! Evaluate the solution's second derivative with respect to \c x
    Scalar _xx(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const;

    //! Evaluate the solution's derivative with respect to \c x and \c y
    Scalar _xy(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const;

    //! Evaluate the solution's derivative with respect to \c x and \c z
    Scalar _xz(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const;

    //! Evaluate the solution's derivative with respect to \c y
    Scalar _y(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const;

    //! Evaluate the solution's second derivative with respect to \c y
    Scalar _yy(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const;

    //! Evaluate the solution's derivative with respect to \c y and \c z
    Scalar _yz(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const;

    //! Evaluate the solution's derivative with respect to \c z
    Scalar _z(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const;

    //! Evaluate the solution's second derivative with respect to \c z
    Scalar _zz(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const;

}; // end class

template<typename Scalar>
Scalar solution<Scalar>::operator()(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const
{
    return a_0*::std::cos(g_0 + f_0*t) + a_x*::std::cos(c_x + b_x*x)*::std::cos(g_x + f_x*t) + a_y*::std::cos(g_y + f_y*t)*::std::cos(c_y + b_y*y) + a_z*::std::cos(g_z + f_z*t)*::std::cos(c_z + b_z*z) + a_xy*::std::cos(g_xy + f_xy*t)*::std::cos(c_xy + b_xy*x)*::std::cos(e_xy + d_xy*y) + a_xz*::std::cos(c_xz + b_xz*x)*::std::cos(g_xz + f_xz*t)*::std::cos(e_xz + d_xz*z) + a_yz*::std::cos(c_yz + b_yz*y)*::std::cos(e_yz + d_yz*z)*::std::cos(g_yz + f_yz*t);
}

template<typename Scalar>
Scalar solution<Scalar>::_t(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const
{
    return -a_0*f_0*::std::sin(g_0 + f_0*t) - a_x*f_x*::std::cos(c_x + b_x*x)*::std::sin(g_x + f_x*t) - a_y*f_y*::std::cos(c_y + b_y*y)*::std::sin(g_y + f_y*t) - a_z*f_z*::std::cos(c_z + b_z*z)*::std::sin(g_z + f_z*t) - a_xy*f_xy*::std::cos(c_xy + b_xy*x)*::std::cos(e_xy + d_xy*y)*::std::sin(g_xy + f_xy*t) - a_xz*f_xz*::std::cos(c_xz + b_xz*x)*::std::cos(e_xz + d_xz*z)*::std::sin(g_xz + f_xz*t) - a_yz*f_yz*::std::cos(c_yz + b_yz*y)*::std::cos(e_yz + d_yz*z)*::std::sin(g_yz + f_yz*t);
}

template<typename Scalar>
Scalar solution<Scalar>::_x(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const
{
    return -a_x*b_x*::std::cos(g_x + f_x*t)*::std::sin(c_x + b_x*x) - a_xy*b_xy*::std::cos(g_xy + f_xy*t)*::std::cos(e_xy + d_xy*y)*::std::sin(c_xy + b_xy*x) - a_xz*b_xz*::std::cos(g_xz + f_xz*t)*::std::cos(e_xz + d_xz*z)*::std::sin(c_xz + b_xz*x);
}

template<typename Scalar>
Scalar solution<Scalar>::_xx(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const
{
    return -a_x*::std::pow(b_x,static_cast<Scalar>(2))*::std::cos(c_x + b_x*x)*::std::cos(g_x + f_x*t) - a_xy*::std::pow(b_xy,static_cast<Scalar>(2))*::std::cos(g_xy + f_xy*t)*::std::cos(c_xy + b_xy*x)*::std::cos(e_xy + d_xy*y) - a_xz*::std::pow(b_xz,static_cast<Scalar>(2))*::std::cos(c_xz + b_xz*x)*::std::cos(g_xz + f_xz*t)*::std::cos(e_xz + d_xz*z);
}

template<typename Scalar>
Scalar solution<Scalar>::_xy(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const
{
    return a_xy*b_xy*d_xy*::std::cos(g_xy + f_xy*t)*::std::sin(c_xy + b_xy*x)*::std::sin(e_xy + d_xy*y);
}

template<typename Scalar>
Scalar solution<Scalar>::_xz(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const
{
    return a_xz*b_xz*d_xz*::std::cos(g_xz + f_xz*t)*::std::sin(c_xz + b_xz*x)*::std::sin(e_xz + d_xz*z);
}

template<typename Scalar>
Scalar solution<Scalar>::_y(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const
{
    return -a_y*b_y*::std::cos(g_y + f_y*t)*::std::sin(c_y + b_y*y) - a_xy*d_xy*::std::cos(g_xy + f_xy*t)*::std::cos(c_xy + b_xy*x)*::std::sin(e_xy + d_xy*y) - a_yz*b_yz*::std::cos(e_yz + d_yz*z)*::std::cos(g_yz + f_yz*t)*::std::sin(c_yz + b_yz*y);
}

template<typename Scalar>
Scalar solution<Scalar>::_yy(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const
{
    return -a_y*::std::pow(b_y,static_cast<Scalar>(2))*::std::cos(g_y + f_y*t)*::std::cos(c_y + b_y*y) - a_xy*::std::pow(d_xy,static_cast<Scalar>(2))*::std::cos(g_xy + f_xy*t)*::std::cos(c_xy + b_xy*x)*::std::cos(e_xy + d_xy*y) - a_yz*::std::pow(b_yz,static_cast<Scalar>(2))*::std::cos(c_yz + b_yz*y)*::std::cos(e_yz + d_yz*z)*::std::cos(g_yz + f_yz*t);
}

template<typename Scalar>
Scalar solution<Scalar>::_yz(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const
{
    return a_yz*b_yz*d_yz*::std::cos(g_yz + f_yz*t)*::std::sin(c_yz + b_yz*y)*::std::sin(e_yz + d_yz*z);
}

template<typename Scalar>
Scalar solution<Scalar>::_z(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const
{
    return -a_z*b_z*::std::cos(g_z + f_z*t)*::std::sin(c_z + b_z*z) - a_xz*d_xz*::std::cos(c_xz + b_xz*x)*::std::cos(g_xz + f_xz*t)*::std::sin(e_xz + d_xz*z) - a_yz*d_yz*::std::cos(c_yz + b_yz*y)*::std::cos(g_yz + f_yz*t)*::std::sin(e_yz + d_yz*z);
}

template<typename Scalar>
Scalar solution<Scalar>::_zz(const Scalar x, const Scalar y, const Scalar z, const Scalar t) const
{
    return -a_z*::std::pow(b_z,static_cast<Scalar>(2))*::std::cos(g_z + f_z*t)*::std::cos(c_z + b_z*z) - a_xz*::std::pow(d_xz,static_cast<Scalar>(2))*::std::cos(c_xz + b_xz*x)*::std::cos(g_xz + f_xz*t)*::std::cos(e_xz + d_xz*z) - a_yz*::std::pow(d_yz,static_cast<Scalar>(2))*::std::cos(c_yz + b_yz*y)*::std::cos(e_yz + d_yz*z)*::std::cos(g_yz + f_yz*t);
}
