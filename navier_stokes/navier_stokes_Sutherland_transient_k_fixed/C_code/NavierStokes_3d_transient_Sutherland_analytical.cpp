double eval_exact_rho(double x, double y, double z, double t)
{
  using std::cos;
  using std::sin;
  
  double rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi * t / Lt);
  
  return rho_an;
}

double eval_exact_p(double x, double y, double z, double t)
{
  using std::cos;
  using std::sin;
  
  double p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt);
  
  return p_an;
}

double eval_exact_u(double x, double y, double z, double t)
{
  using std::cos;
  using std::sin;
  
  double u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt);
  
  return u_an;
}

double eval_exact_v(double x, double y, double z, double t)
{
  using std::cos;
  using std::sin;
  
  double v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt);
  
  return v_an;
}


double eval_exact_w(double x, double y, double z, double t)
{
  using std::cos;
  using std::sin;
  
  double w_an = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt);
  
  return w_an;
}




