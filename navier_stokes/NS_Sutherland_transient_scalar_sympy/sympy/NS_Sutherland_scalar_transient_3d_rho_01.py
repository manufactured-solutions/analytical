# This program calculates the source term Q for the 3D transient Navier-Stokes equations with Sutherland viscosity model: 
# mu=A_mu*T^(3/2)/(T+B_mu).

# Density rho


from sympy import *

var('x y z t')

u = Function('u')(x,y,z,t)
v = Function('v')(x,y,z,t)
w = Function('w')(x,y,z,t)
rho = Function('rho')(x,y,z,t)


#3D Transient Navier-Stokes equation - density rho:
diff(rho, t)+diff(rho*u, x)+diff(rho*v, y)+diff(rho*w, z)==0;


# Mass conservation equation, writen as differential operator:
Lo=diff(rho, t)+diff(rho*u, x)+diff(rho*v, y)+diff(rho*w, z);

L1=diff(rho*u, x)+diff(rho*v, y)+diff(rho*w, z)
L2=diff(rho, t)
 
 
Res=Lo-(L1+L2)
Res=expand(Res)
Res== 0 #True

# Manufactured solutions ----------------------------------------------------------------------------------------
var("""L,Lt,rho_0,rho_x,a_rhox,rho_y,a_rhoy,rho_z,a_rhoz,rho_t,a_rhot, u_0,u_x,a_ux,u_y,a_uy,u_z,a_uz,u_t,a_ut,
v_0,v_x,a_vx,v_y,a_vy,v_z,a_vz,v_t,a_vt,w_0,w_x,a_wx,w_y,a_wy,w_z,a_wz,w_t,a_wt,p_0,p_x,a_px,p_y,a_py,p_z,a_pz,p_t,a_pt,""",real=True)

rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi * t / Lt);
u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt);
v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt);
w_an = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt);
p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt);


# Aplying L1 on rho in order to obtain source term Q1 ----------------------
Q1=L1.subs(u,u_an)
Q1=Q1.subs(v,v_an)
Q1=Q1.subs(w,w_an)
Q1=Q1.subs(rho,rho_an)



# Aplying L2 on rho, u, v and w in order to obtain source term Q2 ----------
Q2=L2.subs(u,u_an)
Q2=Q2.subs(v,v_an)
Q2=Q2.subs(w,w_an)
Q2=Q2.subs(rho,rho_an)



#--------------------------------------------------------------------------
# Source tem for the continuity equation
Qrho=Q1+Q2



