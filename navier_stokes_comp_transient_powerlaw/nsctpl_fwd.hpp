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
 * standards and well-aligned with mathematical notation.  The
 * foreach_parameter() member allows invoking an operation on each solution
 * parameter to aid parameter registration, initialization, or output.
 * Parameters may have an infix name added for use with foreach_parameter().
 * For example, setting the name 'phi' will cause foreach_parameter() to report
 * names like 'a_phix'.
 */
template<typename Scalar>
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
    template<typename BinaryFunction>
    void foreach_parameter(BinaryFunction& f) const {
        ::std::ostringstream os;
        FOR_ALL_SOLUTION_PARAMETERS(APPLY)
    }

    //! Invoke the binary function f on each parameter name and its value.
    template<typename BinaryFunction>
    void foreach_parameter(BinaryFunction& f) {
        ::std::ostringstream os;
        FOR_ALL_SOLUTION_PARAMETERS(APPLY)
    }

#undef APPLY
#undef APPLY_STRINGIFY

// Avoid polluting the global namespace with an implementation detail
#undef FOR_ALL_SOLUTION_PARAMETERS

    //! Evaluate the solution
    Scalar operator()(Scalar x, Scalar y, Scalar z, Scalar t) const;

    //! Evaluate the solution's derivative with respect to time
    Scalar _t(Scalar x, Scalar y, Scalar z, Scalar t) const;

    //! Evaluate the solution's derivative with respect to \c x
    Scalar _x(Scalar x, Scalar y, Scalar z, Scalar t) const;

    //! Evaluate the solution's second derivative with respect to \c x
    Scalar _xx(Scalar x, Scalar y, Scalar z, Scalar t) const;

    //! Evaluate the solution's derivative with respect to \c x and \c y
    Scalar _xy(Scalar x, Scalar y, Scalar z, Scalar t) const;

    //! Evaluate the solution's derivative with respect to \c x and \c z
    Scalar _xz(Scalar x, Scalar y, Scalar z, Scalar t) const;

    //! Evaluate the solution's derivative with respect to \c y
    Scalar _y(Scalar x, Scalar y, Scalar z, Scalar t) const;

    //! Evaluate the solution's second derivative with respect to \c y
    Scalar _yy(Scalar x, Scalar y, Scalar z, Scalar t) const;

    //! Evaluate the solution's derivative with respect to \c y and \c z
    Scalar _yz(Scalar x, Scalar y, Scalar z, Scalar t) const;

    //! Evaluate the solution's derivative with respect to \c z
    Scalar _z(Scalar x, Scalar y, Scalar z, Scalar t) const;

    //! Evaluate the solution's second derivative with respect to \c z
    Scalar _zz(Scalar x, Scalar y, Scalar z, Scalar t) const;

}; // end class


/**
 * Template that, given a primitive solution for \c rho, \c u, \c v, \c w, and
 * \c T along with a floating point type, evaluates a manufactured solution and
 * the required forcing for the transient, compressible Navier--Stokes
 * equations with a power law viscosity.
 */
template<template <typename T> class PrimitiveSolution, typename Scalar>
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
    PrimitiveSolution<Scalar> soln_rho;  //!< Analytic solution for rho
    PrimitiveSolution<Scalar> soln_u;    //!< Analytic solution for u
    PrimitiveSolution<Scalar> soln_v;    //!< Analytic solution for v
    PrimitiveSolution<Scalar> soln_w;    //!< Analytic solution for w
    PrimitiveSolution<Scalar> soln_T;    //!< Analytic solution for T
    //!@}

    //! Default constructor
    generic_manufactured_solution()
        : gamma(0), R(0), beta(0), mu_r(0), T_r(0), k_r(0), lambda_r(0),
          soln_rho("rho"), soln_u("u"), soln_v("v"), soln_w("w"), soln_T("T")
    {
    }

    //! Invoke the binary function f on each parameter name and its value.
    template<typename BinaryFunction>
    void foreach_parameter(BinaryFunction& f) const {
        f(::std::string("gamma"),    gamma   );
        f(::std::string("R"),        R       );
        f(::std::string("beta"),     beta    );
        f(::std::string("mu_r"),     mu_r    );
        f(::std::string("T_r"),      T_r     );
        f(::std::string("k_r"),      k_r     );
        f(::std::string("lambda_r"), lambda_r);
        soln_rho.foreach_parameter(f);
        soln_u.foreach_parameter(f);
        soln_v.foreach_parameter(f);
        soln_w.foreach_parameter(f);
        soln_T.foreach_parameter(f);
    }

    //! Invoke the binary function f on each parameter name and its value.
    template<typename BinaryFunction>
    void foreach_parameter(BinaryFunction& f) {
        f(::std::string("gamma"),    gamma   );
        f(::std::string("R"),        R       );
        f(::std::string("beta"),     beta    );
        f(::std::string("mu_r"),     mu_r    );
        f(::std::string("T_r"),      T_r     );
        f(::std::string("k_r"),      k_r     );
        f(::std::string("lambda_r"), lambda_r);
        soln_rho.foreach_parameter(f);
        soln_u.foreach_parameter(f);
        soln_v.foreach_parameter(f);
        soln_w.foreach_parameter(f);
        soln_T.foreach_parameter(f);
    }

    // Analytically determined quantities
    Scalar eval_exact_rho(Scalar x, Scalar y, Scalar z, Scalar t) const;
    Scalar eval_exact_u  (Scalar x, Scalar y, Scalar z, Scalar t) const;
    Scalar eval_exact_v  (Scalar x, Scalar y, Scalar z, Scalar t) const;
    Scalar eval_exact_w  (Scalar x, Scalar y, Scalar z, Scalar t) const;
    Scalar eval_exact_T  (Scalar x, Scalar y, Scalar z, Scalar t) const;
    Scalar eval_g_rho    (Scalar x, Scalar y, Scalar z, Scalar t, int direction) const;
    Scalar eval_g_u      (Scalar x, Scalar y, Scalar z, Scalar t, int direction) const;
    Scalar eval_g_v      (Scalar x, Scalar y, Scalar z, Scalar t, int direction) const;
    Scalar eval_g_w      (Scalar x, Scalar y, Scalar z, Scalar t, int direction) const;
    Scalar eval_g_T      (Scalar x, Scalar y, Scalar z, Scalar t, int direction) const;

    // Quantities built from the analytical solutions
    // TODO Build up eval_q_u, eval_q_v, eval_q_w, eval_q_e, eval_q_T, eval_q_p
    Scalar eval_exact_e  (Scalar x, Scalar y, Scalar z, Scalar t) const;
    Scalar eval_exact_p  (Scalar x, Scalar y, Scalar z, Scalar t) const;
    Scalar eval_exact_mu (Scalar x, Scalar y, Scalar z, Scalar t) const;
    Scalar eval_g_e      (Scalar x, Scalar y, Scalar z, Scalar t, int direction) const;
    Scalar eval_g_p      (Scalar x, Scalar y, Scalar z, Scalar t, int direction) const;
    Scalar eval_g_mu     (Scalar x, Scalar y, Scalar z, Scalar t, int direction) const;
    Scalar eval_q_rho    (Scalar x, Scalar y, Scalar z, Scalar t) const;
    Scalar eval_q_rho_u  (Scalar x, Scalar y, Scalar z, Scalar t) const;
    Scalar eval_q_rho_v  (Scalar x, Scalar y, Scalar z, Scalar t) const;
    Scalar eval_q_rho_w  (Scalar x, Scalar y, Scalar z, Scalar t) const;
    Scalar eval_q_rho_e  (Scalar x, Scalar y, Scalar z, Scalar t) const;

}; // end class

/**
 * Template that, given a floating point type, evaluates a manufactured
 * solution and the required forcing for the transient, compressible
 * Navier--Stokes equations with a power law viscosity.
 */
template<typename Scalar>
class manufactured_solution
    : public generic_manufactured_solution<primitive_solution,Scalar> {};

} // end namespace nsctpl
