u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L);
v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L);
u_an_t = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_t * cos(a_ut * pi * t / L);
v_an_t = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_t * sin(a_vt * pi * t / L);
grad_u_an[0] = u_x * cos(a_ux * pi * x / L) * a_ux * pi / L;
grad_u_an[1] = -u_y * sin(a_uy * pi * y / L) * a_uy * pi / L;
grad_v_an[0] = -v_x * sin(a_vx * pi * x / L) * a_vx * pi / L;
grad_v_an[1] = v_y * cos(a_vy * pi * y / L) * a_vy * pi / L;
