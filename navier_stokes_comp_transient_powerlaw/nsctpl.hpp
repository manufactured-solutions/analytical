template<typename Scalar>
Scalar nsctpl_solution<Scalar>::operator()(Scalar x, Scalar y, Scalar z, Scalar t) const
{
    return a_0*::std::cos(g_0 + f_0*t) + a_x*::std::cos(c_x + b_x*x)*::std::cos(g_x + f_x*t) + a_y*::std::cos(g_y + f_y*t)*::std::cos(c_y + b_y*y) + a_z*::std::cos(g_z + f_z*t)*::std::cos(c_z + b_z*z) + a_xy*::std::cos(g_xy + f_xy*t)*::std::cos(c_xy + b_xy*x)*::std::cos(e_xy + d_xy*y) + a_xz*::std::cos(c_xz + b_xz*x)*::std::cos(g_xz + f_xz*t)*::std::cos(e_xz + d_xz*z) + a_yz*::std::cos(c_yz + b_yz*y)*::std::cos(e_yz + d_yz*z)*::std::cos(g_yz + f_yz*t);
}

template<typename Scalar>
Scalar nsctpl_solution<Scalar>::_t(Scalar x, Scalar y, Scalar z, Scalar t) const
{
    return -a_0*f_0*::std::sin(g_0 + f_0*t) - a_x*f_x*::std::cos(c_x + b_x*x)*::std::sin(g_x + f_x*t) - a_y*f_y*::std::cos(c_y + b_y*y)*::std::sin(g_y + f_y*t) - a_z*f_z*::std::cos(c_z + b_z*z)*::std::sin(g_z + f_z*t) - a_xy*f_xy*::std::cos(c_xy + b_xy*x)*::std::cos(e_xy + d_xy*y)*::std::sin(g_xy + f_xy*t) - a_xz*f_xz*::std::cos(c_xz + b_xz*x)*::std::cos(e_xz + d_xz*z)*::std::sin(g_xz + f_xz*t) - a_yz*f_yz*::std::cos(c_yz + b_yz*y)*::std::cos(e_yz + d_yz*z)*::std::sin(g_yz + f_yz*t);
}

template<typename Scalar>
Scalar nsctpl_solution<Scalar>::_x(Scalar x, Scalar y, Scalar z, Scalar t) const
{
    return -a_x*b_x*::std::cos(g_x + f_x*t)*::std::sin(c_x + b_x*x) - a_xy*b_xy*::std::cos(g_xy + f_xy*t)*::std::cos(e_xy + d_xy*y)*::std::sin(c_xy + b_xy*x) - a_xz*b_xz*::std::cos(g_xz + f_xz*t)*::std::cos(e_xz + d_xz*z)*::std::sin(c_xz + b_xz*x);
}

template<typename Scalar>
Scalar nsctpl_solution<Scalar>::_xx(Scalar x, Scalar y, Scalar z, Scalar t) const
{
    return -a_x*::std::pow(b_x,static_cast<Scalar>(2))*::std::cos(c_x + b_x*x)*::std::cos(g_x + f_x*t) - a_xy*::std::pow(b_xy,static_cast<Scalar>(2))*::std::cos(g_xy + f_xy*t)*::std::cos(c_xy + b_xy*x)*::std::cos(e_xy + d_xy*y) - a_xz*::std::pow(b_xz,static_cast<Scalar>(2))*::std::cos(c_xz + b_xz*x)*::std::cos(g_xz + f_xz*t)*::std::cos(e_xz + d_xz*z);
}

template<typename Scalar>
Scalar nsctpl_solution<Scalar>::_xy(Scalar x, Scalar y, Scalar z, Scalar t) const
{
    return a_xy*b_xy*d_xy*::std::cos(g_xy + f_xy*t)*::std::sin(c_xy + b_xy*x)*::std::sin(e_xy + d_xy*y);
}

template<typename Scalar>
Scalar nsctpl_solution<Scalar>::_xz(Scalar x, Scalar y, Scalar z, Scalar t) const
{
    return a_xz*b_xz*d_xz*::std::cos(g_xz + f_xz*t)*::std::sin(c_xz + b_xz*x)*::std::sin(e_xz + d_xz*z);
}

template<typename Scalar>
Scalar nsctpl_solution<Scalar>::_y(Scalar x, Scalar y, Scalar z, Scalar t) const
{
    return -a_y*b_y*::std::cos(g_y + f_y*t)*::std::sin(c_y + b_y*y) - a_xy*d_xy*::std::cos(g_xy + f_xy*t)*::std::cos(c_xy + b_xy*x)*::std::sin(e_xy + d_xy*y) - a_yz*b_yz*::std::cos(e_yz + d_yz*z)*::std::cos(g_yz + f_yz*t)*::std::sin(c_yz + b_yz*y);
}

template<typename Scalar>
Scalar nsctpl_solution<Scalar>::_yy(Scalar x, Scalar y, Scalar z, Scalar t) const
{
    return -a_y*::std::pow(b_y,static_cast<Scalar>(2))*::std::cos(g_y + f_y*t)*::std::cos(c_y + b_y*y) - a_xy*::std::pow(d_xy,static_cast<Scalar>(2))*::std::cos(g_xy + f_xy*t)*::std::cos(c_xy + b_xy*x)*::std::cos(e_xy + d_xy*y) - a_yz*::std::pow(b_yz,static_cast<Scalar>(2))*::std::cos(c_yz + b_yz*y)*::std::cos(e_yz + d_yz*z)*::std::cos(g_yz + f_yz*t);
}

template<typename Scalar>
Scalar nsctpl_solution<Scalar>::_yz(Scalar x, Scalar y, Scalar z, Scalar t) const
{
    return a_yz*b_yz*d_yz*::std::cos(g_yz + f_yz*t)*::std::sin(c_yz + b_yz*y)*::std::sin(e_yz + d_yz*z);
}

template<typename Scalar>
Scalar nsctpl_solution<Scalar>::_z(Scalar x, Scalar y, Scalar z, Scalar t) const
{
    return -a_z*b_z*::std::cos(g_z + f_z*t)*::std::sin(c_z + b_z*z) - a_xz*d_xz*::std::cos(c_xz + b_xz*x)*::std::cos(g_xz + f_xz*t)*::std::sin(e_xz + d_xz*z) - a_yz*d_yz*::std::cos(c_yz + b_yz*y)*::std::cos(g_yz + f_yz*t)*::std::sin(e_yz + d_yz*z);
}

template<typename Scalar>
Scalar nsctpl_solution<Scalar>::_zz(Scalar x, Scalar y, Scalar z, Scalar t) const
{
    return -a_z*::std::pow(b_z,static_cast<Scalar>(2))*::std::cos(g_z + f_z*t)*::std::cos(c_z + b_z*z) - a_xz*::std::pow(d_xz,static_cast<Scalar>(2))*::std::cos(c_xz + b_xz*x)*::std::cos(g_xz + f_xz*t)*::std::cos(e_xz + d_xz*z) - a_yz*::std::pow(d_yz,static_cast<Scalar>(2))*::std::cos(c_yz + b_yz*y)*::std::cos(e_yz + d_yz*z)*::std::cos(g_yz + f_yz*t);
}
