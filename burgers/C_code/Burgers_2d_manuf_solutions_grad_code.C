u_an = u_0 * (sin(x * x + y * y + omega * t) + epsilon);
v_an = v_0 * (cos(x * x + y * y + omega * t) + epsilon);
grad_u_an[0] = 0.2e1 * u_0 * cos(x * x + y * y + omega * t) * x;
grad_u_an[1] = 0.2e1 * u_0 * cos(x * x + y * y + omega * t) * y;
grad_v_an[0] = -0.2e1 * v_0 * sin(x * x + y * y + omega * t) * x;
grad_v_an[1] = -0.2e1 * v_0 * sin(x * x + y * y + omega * t) * y;
