// Manufactured solutions for the 2D steady Burgers equations and their gradients
void Steady_Manuf_Solutions (double x, double y)
{
  double u_an;
  double v_an;
  double grad_u_an[2];
  double grad_v_an[2];
    // "Manufactured solutions, steady -------------------------------------"
;
  u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
  v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
    // "Gradients of Manufactured solutions -------------------------"
;
  grad_u_an[0] = u_x * cos(a_ux * PI * x / L) * a_ux * PI / L;
  grad_u_an[1] = -u_y * sin(a_uy * PI * y / L) * a_uy * PI / L;
  grad_v_an[0] = -v_x * sin(a_vx * PI * x / L) * a_vx * PI / L;
  grad_v_an[1] = v_y * cos(a_vy * PI * y / L) * a_vy * PI / L;
}
// Manufactured solutions for the 2D transient Burgers equations and their gradients
void Transient_Manuf_Solutions (double x, double y, double t)
{
  double u_an;
  double v_an;
  double grad_u_an[2];
  double grad_v_an[2];
    // "Manufactured solutions, transient -------------------------------------"
;
  u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_t * cos(a_ut * pi * t / L);
  v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_t * sin(a_vt * pi * t / L);
    // "Gradients of Manufactured solutions -------------------------"
;
  grad_u_an[0] = u_x * cos(a_ux * PI * x / L) * a_ux * PI / L;
  grad_u_an[1] = -u_y * sin(a_uy * PI * y / L) * a_uy * PI / L;
  grad_v_an[0] = -v_x * sin(a_vx * PI * x / L) * a_vx * PI / L;
  grad_v_an[1] = v_y * cos(a_vy * PI * y / L) * a_vy * PI / L;
}
