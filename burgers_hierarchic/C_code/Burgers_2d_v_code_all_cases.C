// Source term for the y-momentum equation (v), 2 D transient, viscous Burgers equations
double SourceQ_v_transient_viscous (
  double x,
  double y,
  double t,
  double nu)
{
  double Qv_tv;
  double U;
  double V;
  double Q_v_time;
  double Q_v_convection;
  double Q_v_dissipation;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / Lt);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / Lt);
    // "Contribution from the time derivative  to the total source term ----------------------------------------------"
;
  Q_v_time = a_vt * PI * v_t * cos(a_vt * PI * t / Lt) / Lt;
    // "Contribution from the convective terms to the total source term ---------------------------------------------"
;
  Q_v_convection = -a_vx * PI * v_x * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * V / L;
    // "Contribution from the viscous/dissipation terms to the total source term ----------------------------------------"
;
  Q_v_dissipation = a_vx * a_vx * PI * PI * v_x * nu * cos(a_vx * PI * x / L) * pow(L, -0.2e1) + a_vy * a_vy * PI * PI * v_y * nu * sin(a_vy * PI * y / L) * pow(L, -0.2e1);
    // "Total source term ------------------------------------------------------------------------------------------"
;
  Qv_tv = Q_v_dissipation + Q_v_convection + Q_v_time;
  return(Qv_tv);
}
// Source term for the y-momentum equation (v), 2 D transient, viscous Burgers equations
double SourceQ_v_transient_inviscid (double x, double y, double t)
{
  double Qv_tinv;
  double U;
  double V;
  double Q_v_time;
  double Q_v_convection;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L) + u_t * cos(a_ut * PI * t / Lt);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L) + v_t * sin(a_vt * PI * t / Lt);
    // "Contribution from the time derivative  to the total source term ----------------------------------------------"
;
  Q_v_time = a_vt * PI * v_t * cos(a_vt * PI * t / Lt) / Lt;
    // "Contribution from the convective terms to the total source term ---------------------------------------------"
;
  Q_v_convection = -a_vx * PI * v_x * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * V / L;
    // "Total source term ------------------------------------------------------------------------------------------"
;
  Qv_tinv = Q_v_convection + Q_v_time;
  return(Qv_tinv);
}
// Source term for the y-momentum equation (v), 2D steady, viscous Burgers equations
double SourceQ_v_steady_viscous (double x, double y, double nu)
{
  double Qv_sv;
  double U;
  double V;
  double Q_v_convection;
  double Q_v_dissipation;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);
    // "Contribution from the convective terms to the total source term ---------------------------------------------"
;
  Q_v_convection = -a_vx * PI * v_x * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * V / L;
    // "Contribution from the viscous/dissipation terms to the total source term ----------------------------------------"
;
  Q_v_dissipation = a_vx * a_vx * PI * PI * v_x * nu * cos(a_vx * PI * x / L) * pow(L, -0.2e1) + a_vy * a_vy * PI * PI * v_y * nu * sin(a_vy * PI * y / L) * pow(L, -0.2e1);
    // "Total source term ------------------------------------------------------------------------------------------"
;
  Qv_sv = Q_v_dissipation + Q_v_convection;
  return(Qv_sv);
}
// Source term for the y-momentum equation (v), 2D steady, inviscid Burgers equations
double SourceQ_v_steady_inviscid (double x, double y)
{
  double Qv_sinv;
  double U;
  double V;
  double Q_v_convection;
  U = u_0 + u_x * sin(a_ux * PI * x / L) + u_y * cos(a_uy * PI * y / L);
  V = v_0 + v_x * cos(a_vx * PI * x / L) + v_y * sin(a_vy * PI * y / L);
    // "Contribution from the convective terms to the total source term ---------------------------------------------"
;
  Q_v_convection = -a_vx * PI * v_x * U * sin(a_vx * PI * x / L) / L + (a_ux * u_x * cos(a_ux * PI * x / L) + 0.2e1 * a_vy * v_y * cos(a_vy * PI * y / L)) * PI * V / L;
    // "Total source term ------------------------------------------------------------------------------------------"
;
  Qv_sinv = Q_v_convection;
  return(Qv_sinv);
}
