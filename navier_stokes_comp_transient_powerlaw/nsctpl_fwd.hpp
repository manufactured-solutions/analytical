// Using this file will require #include-ing <cmath>, <limits>, <sstream>, and
// <string>.  These headers are not included here to allow #include-ing this
// file inside any namespace, as is the consistent use of ::std instead of
// merely std.

namespace nsctpl {

/**
 * Class providing operations for evaluating an analytical solution at a
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
 * standards and well-aligned with mathematical notation.  Location and
 * time arguments are templated to allow extended precision intermediate
 * computations (followed by truncation) when possible.
 *
 * The foreach_parameter() member allows invoking an operation on each solution
 * parameter to aid parameter registration, initialization, or output.
 * Parameters may have an infix name added for use with foreach_parameter().
 * For example, setting the name 'phi' will cause foreach_parameter() to report
 * names like 'a_phix'.
 */
template <typename Scalar>
class primitive_solution {

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

    //! The name used infix in foreach_parameter' names.  For example,
    //! name = 'phi' implies parameters names like 'a_phix'.
    const ::std::string name;

    //! Construct an instance using \c name in the reported parameter names.
    //! All parameters set to zero at construction time.
    explicit primitive_solution(const ::std::string &name = "")
#define APPLY(pre,suf) pre##suf(0),
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
    template <typename BinaryFunction>
    void foreach_parameter(BinaryFunction& f) const {
        ::std::ostringstream os;
        FOR_ALL_SOLUTION_PARAMETERS(APPLY)
    }

    //! Invoke the binary function f on each parameter name and its value.
    template <typename BinaryFunction>
    void foreach_parameter(BinaryFunction& f) {
        ::std::ostringstream os;
        FOR_ALL_SOLUTION_PARAMETERS(APPLY)
    }

#undef APPLY
#undef APPLY_STRINGIFY

// Avoid polluting the global namespace with an implementation detail
#undef FOR_ALL_SOLUTION_PARAMETERS

    //! Evaluate the solution
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar operator()(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to time
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _t(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c x
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _x(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's second derivative with respect to \c x
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _xx(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c x and \c y
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _xy(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c x and \c z
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _xz(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c y
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _y(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's second derivative with respect to \c y
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _yy(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c y and \c z
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _yz(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's derivative with respect to \c z
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _z(T1 x, T2 y, T3 z, T4 t) const;

    //! Evaluate the solution's second derivative with respect to \c z
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar _zz(T1 x, T2 y, T3 z, T4 t) const;

}; // end class


/**
 * Template that, given a primitive solution for \c rho, \c u, \c v, \c w, and
 * \c T along with a floating point type, evaluates a manufactured solution and
 * the required forcing for the transient, compressible Navier--Stokes
 * equations with a power law viscosity.  Location and time arguments are
 * templated to allow extended precision intermediate computations (followed by
 * truncation) when possible.
 */
template <template <typename T> class PrimitiveSolution, typename Scalar>
class generic_manufactured_solution {

public:

    //! Scenario parameters
    //!@{
    Scalar gamma;                    //!< Constant ratio of specific heats
    Scalar R;                        //!< Gas constant
    Scalar beta;                     //!< Viscosity power law exponent
    Scalar mu_r;                     //!< Reference dynamic viscosity
    Scalar T_r;                      //!< Reference temperature
    Scalar k_r;                      //!< Reference thermal conductivity
    Scalar lambda_r;                 //!< Reference bulk viscosity
    //!@}

    //! Analytic solutions (which contain additional parameters)
    //!@{
    PrimitiveSolution<Scalar> rho;  //!< Analytic solution for rho
    PrimitiveSolution<Scalar> u;    //!< Analytic solution for u
    PrimitiveSolution<Scalar> v;    //!< Analytic solution for v
    PrimitiveSolution<Scalar> w;    //!< Analytic solution for w
    PrimitiveSolution<Scalar> T;    //!< Analytic solution for T
    //!@}

    //! Default constructor
    generic_manufactured_solution()
        : gamma(0), R(0), beta(0), mu_r(0), T_r(0), k_r(0), lambda_r(0),
          rho("rho"), u("u"), v("v"), w("w"), T("T")
    {
    }

    //! Invoke the binary function f on each parameter name and its value.
    template <typename BinaryFunction>
    void foreach_parameter(BinaryFunction& f) const {
        f(::std::string("gamma"),    gamma   );
        f(::std::string("R"),        R       );
        f(::std::string("beta"),     beta    );
        f(::std::string("mu_r"),     mu_r    );
        f(::std::string("T_r"),      T_r     );
        f(::std::string("k_r"),      k_r     );
        f(::std::string("lambda_r"), lambda_r);
        rho.foreach_parameter(f);
        u.foreach_parameter(f);
        v.foreach_parameter(f);
        w.foreach_parameter(f);
        T.foreach_parameter(f);
    }

    //! Invoke the binary function f on each parameter name and its value.
    template <typename BinaryFunction>
    void foreach_parameter(BinaryFunction& f) {
        f(::std::string("gamma"),    gamma   );
        f(::std::string("R"),        R       );
        f(::std::string("beta"),     beta    );
        f(::std::string("mu_r"),     mu_r    );
        f(::std::string("T_r"),      T_r     );
        f(::std::string("k_r"),      k_r     );
        f(::std::string("lambda_r"), lambda_r);
        rho.foreach_parameter(f);
        u.foreach_parameter(f);
        v.foreach_parameter(f);
        w.foreach_parameter(f);
        T.foreach_parameter(f);
    }

    // Analytically determined quantities
    // Note that PrimitiveSolution members can be used directly.
    // For example, T(x,y,z,t) or T._xx(x,y,z,t)
    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_rho (T1 x, T2 y, T3 z, T4 t, int direction) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_u   (T1 x, T2 y, T3 z, T4 t, int direction) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_v   (T1 x, T2 y, T3 z, T4 t, int direction) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_w   (T1 x, T2 y, T3 z, T4 t, int direction) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_T   (T1 x, T2 y, T3 z, T4 t, int direction) const;

    // Quantities built from the analytical solutions
    // TODO Build up q_u, q_v, q_w, q_e, q_T, q_p

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar e(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar p(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar mu(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar rhou(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar rhov(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar rhow(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar rhoe(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_e(T1 x, T2 y, T3 z, T4 t, int direction) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_p(T1 x, T2 y, T3 z, T4 t, int direction) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar grad_mu(T1 x, T2 y, T3 z, T4 t, int direction) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar Q_rho(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar Q_rhou(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar Q_rhov(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar Q_rhow(T1 x, T2 y, T3 z, T4 t) const;

    template <typename T1, typename T2, typename T3, typename T4>
    Scalar Q_rhoe(T1 x, T2 y, T3 z, T4 t) const;

}; // end class

/**
 * Template that, given a floating point type, evaluates a manufactured
 * solution and the required forcing for the transient, compressible
 * Navier--Stokes equations with a power law viscosity.
 */
template <typename Scalar>
class manufactured_solution
    : public generic_manufactured_solution<primitive_solution,Scalar> {};

} // end namespace nsctpl
