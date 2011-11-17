// Source term for the x-momentum equation (u), 2 D transient, viscous Burgers equations
double SourceQ_u_transient_viscous (
  double x,
  double y,
  double t,
  double nu)
{
  double Qu_tv;
  double U;
  double V;
  double Q_u_time;
  double Q_u_convection;
  double Q_u_dissipation;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / Lt);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / Lt);
    // "Contribution from the time derivative  to the total source term ----------------------------------------------"
;
  Q_u_time = -a_ut * PI * u_t * sin(a_ut * PI * t / Lt) / Lt;
    // "Contribution from the convective terms to the total source term ---------------------------------------------"
;
  Q_u_convection = -a_uy * PI * u_y * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * U / L;
    // "Contribution from the viscous/dissipation terms to the total source term ----------------------------------------"
;
  Q_u_dissipation = a_ux * a_ux * PI * PI * u_x * nu * sin(a_ux * PI * x / L) * pow(L, -0.2e1) + a_uy * a_uy * PI * PI * u_y * nu * cos(a_uy * PI * y / L) * pow(L, -0.2e1);
    // "Total source term ------------------------------------------------------------------------------------------"
;
  Qu_tv = Q_u_dissipation + Q_u_convection + Q_u_time;
  return(Qu_tv);
}
// Source term for the x-momentum equation (u), 2D transient, inviscid Burgers equations
double SourceQ_u_transient_inviscid (double x, double y, double t)
{
  double Qu_tinv;
  double U;
  double V;
  double Q_u_time;
  double Q_u_convection;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / Lt);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / Lt);
    // "Contribution from the time derivative  to the total source term ----------------------------------------------"
;
  Q_u_time = -a_ut * PI * u_t * sin(a_ut * PI * t / Lt) / Lt;
    // "Contribution from the convective terms to the total source term ---------------------------------------------"
;
  Q_u_convection = -a_uy * PI * u_y * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * U / L;
    // "Total source term ------------------------------------------------------------------------------------------"
;
  Qu_tinv = Q_u_convection + Q_u_time;
  return(Qu_tinv);
}
// Source term for the x-momentum equation (u), 2D steady viscous Burgers equations
double SourceQ_u_steady_viscous (double x, double y, double nu)
{
  double Qu_sv;
  double U;
  double V;
  double Q_u_convection;
  double Q_u_dissipation;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);
    // "Contribution from the convective terms to the total source term ---------------------------------------------"
;
  Q_u_convection = -a_uy * PI * u_y * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * U / L;
    // "Contribution from the viscous/dissipation terms to the total source term ----------------------------------------"
;
  Q_u_dissipation = a_ux * a_ux * PI * PI * u_x * nu * sin(a_ux * PI * x / L) * pow(L, -0.2e1) + a_uy * a_uy * PI * PI * u_y * nu * cos(a_uy * PI * y / L) * pow(L, -0.2e1);
    // "Total source term ------------------------------------------------------------------------------------------"
;
  Qu_sv = Q_u_dissipation + Q_u_convection;
  return(Qu_sv);
}
// Source term for the x-momentum equation (u), 2D steady, inviscid Burgers equations
double SourceQ_u_steady_inviscid (double x, double y)
{
  double Qu_sinv;
  double U;
  double V;
  double Q_u_convection;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);
    // "Contribution from the convective terms to the total source term ---------------------------------------------"
;
  Q_u_convection = -a_uy * PI * u_y * V * sin(a_uy * PI * y / L) / L + (0.2e1 * a_ux * u_x * cos(a_ux * PI * x / L) + a_vy * v_y * cos(a_vy * PI * y / L)) * PI * U / L;
    // "Total source term ------------------------------------------------------------------------------------------"
;
  Qu_sinv = Q_u_convection;
  return(Qu_sinv);
}
