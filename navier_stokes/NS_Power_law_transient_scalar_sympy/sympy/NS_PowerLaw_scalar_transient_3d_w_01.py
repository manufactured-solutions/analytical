# This program calculates the source term Q for the 3D transient Navier-Stokes equations with Power Law viscosity model: 
# mu = mu_r * (T/T_r)^Beta

# Velocity w (momentum)

from sympy import *

var('x y z t')

u = Function('u')(x,y,z,t)
v = Function('v')(x,y,z,t)
w = Function('w')(x,y,z,t)
rho = Function('rho')(x,y,z,t)
p = Function('p')(x,y,z,t)
T = Function('T')(x,y,z,t)
e_t = Function('e_t')(x,y,z,t)

mu = Function('mu')(x,y,z,t)
tau_xx = Function('tau_xx')(x,y,z,t)
tau_xy = Function('tau_xy')(x,y,z,t)
tau_xz = Function('tau_xz')(x,y,z,t)
tau_yx = Function('tau_yx')(x,y,z,t)
tau_yy = Function('tau_yy')(x,y,z,t)
tau_yz = Function('tau_yz')(x,y,z,t)
tau_zx = Function('tau_zx')(x,y,z,t)
tau_zy = Function('tau_zy')(x,y,z,t)
tau_zz = Function('tau_zz')(x,y,z,t)


#3D transient Navier-Stokes equation - velocity w:
diff(rho*w, t)+diff(rho*u*w, x)+diff(rho*v*w, y)+diff(rho*w*w, z)+diff(p, z)==(diff(tau_zx, x))+(diff(tau_zy, y))+(diff(tau_zz, z))
  



# Auxiliary relations -------------------------------------------------------------------------------------------
var("R,k,gamma,E_t,alpha,mu_ref,T_ref,beta")
T = p/(rho*R);
e = R*T/(gamma-1);
e_t = e+(u*u+v*v+w*w)/2;
#mu = mu_ref * (T/T_ref)^beta;
Lambda = (alpha- Rational(2,3) )*mu;

tau_xx = 2 * mu * diff(u, x) + Lambda * ( diff(u, x) + diff(v, y) + diff(w, z) );
tau_yy = 2 * mu * diff(v, y) + Lambda * ( diff(u, x) + diff(v, y) + diff(w, z) );
tau_zz = 2 * mu * diff(w, z) + Lambda * ( diff(u, x) + diff(v, y) + diff(w, z) );

tau_xy = mu*(diff(u, y)+diff(v, x));
tau_xz = mu*(diff(u, z)+diff(w, x));
tau_yz = mu*(diff(w, y)+diff(v, z));
tau_yx = tau_xy;
tau_zx = tau_xz;
tau_zy = tau_yz;

# Velocity-w equation, writen as differential operator:
Lo = diff(rho*w, t)+diff(rho*u*w, x)+diff(rho*v*w, y)+diff(rho*w*w, z)+diff(p, z)-(diff(tau_zx, x)+diff(tau_zy, y)+diff(tau_zz, z))

L1 = +diff(rho*u*w, x)+diff(rho*v*w, y)+diff(rho*w*w, z)
L2 = +diff(p, z)
L3 = -(diff(tau_zx, x)+diff(tau_zy, y)+diff(tau_zz, z))
L4 = diff(rho*w, t)

Res=Lo-L1-L2-L3-L4;
Res=expand(Res)
Res==0 #True


# Manufactured solutions ----------------------------------------------------------------------------------------
var("""L,Lt,rho_0,rho_x,a_rhox,rho_y,a_rhoy,rho_z,a_rhoz,rho_t,a_rhot, u_0,u_x,a_ux,u_y,a_uy,u_z,a_uz,u_t,a_ut,
v_0,v_x,a_vx,v_y,a_vy,v_z,a_vz,v_t,a_vt,w_0,w_x,a_wx,w_y,a_wy,w_z,a_wz,w_t,a_wt,p_0,p_x,a_px,p_y,a_py,p_z,a_pz,p_t,a_pt,""",real=True)

rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi * t
/ Lt);
u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt);
v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt);
w_an = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt);
p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt);


# Aplying L1 on v in order to obtain source term Q1 ----------
Q1=L1.subs(u,u_an)
Q1=Q1.subs(v,v_an)
Q1=Q1.subs(w,w_an)
Q1=Q1.subs(rho,rho_an)
Q1=Q1.subs(p,p_an)




# Aplying L2 on u, v and w in order to obtain source term Q2 ----------
Q2=L2.subs(u,u_an)
Q2=Q2.subs(v,v_an)
Q2=Q2.subs(w,w_an)
Q2=Q2.subs(rho,rho_an)
Q2=Q2.subs(p,p_an)

# Aplying L3 on u, v and w in order to obtain source term Q3 ----------
Q3=L3.subs(u,u_an)
Q3=Q3.subs(v,v_an)
Q3=Q3.subs(w,w_an)
Q3=Q3.subs(rho,rho_an)
Q3=Q3.subs(p,p_an)


# Aplying L4 on u, v and w in order to obtain source term Q4 ----------
Q4=L4.subs(u,u_an)
Q4=Q4.subs(v,v_an)
Q4=Q4.subs(w,w_an)
Q4=Q4.subs(rho,rho_an)
Q4=Q4.subs(p,p_an)




#--------------------------------------------------------------------------
# Source tem for the velocity-w equation
Qw=Q1+Q2+Q3+Q4


