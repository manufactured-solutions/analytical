# This program calculates the source term Q for the 3D transient Navier-Stokes equations with Power Law viscosity model: 
# mu = mu_r * (T/T_r)^Beta

# Scalar phi


from sympy import *

var('x y z t')
var('Gamma')

u = Function('u')(x,y,z,t)
v = Function('v')(x,y,z,t)
w = Function('w')(x,y,z,t)
rho = Function('rho')(x,y,z,t)
phi = Function('phi')(x,y,z,t)

#Passive scalar transport with the 3D Transient Navier-Stokes equation:
diff(rho*phi, t) + diff(rho*u*phi, x)+diff(rho*v*phi, y)+diff(rho*w*phi, z)== Gamma*( diff(phi,x,2)+ diff(phi,y,2) + diff(phi,z,2));


# Mass conservation equation, writen as differential operator:
Lo=diff(rho*phi, t) + diff(rho*u*phi, x)+diff(rho*v*phi, y)+diff(rho*w*phi, z)- Gamma*( diff(phi,x,2)+ diff(phi,y,2) + diff(phi,z,2));

L1=+ diff(rho*u*phi, x)+diff(rho*v*phi, y)+diff(rho*w*phi, z)
L2=diff(rho*phi, t)
L3=- Gamma*( diff(phi,x,2)+ diff(phi,y,2) + diff(phi,z,2)) 

Res=Lo-(L1+L2+L3)
Res=expand(Res)
Res== 0 #True

# Manufactured solutions ----------------------------------------------------------------------------------------
var("""L,Lt,rho_0,rho_x,a_rhox,rho_y,a_rhoy,rho_z,a_rhoz,rho_t,a_rhot, u_0,u_x,a_ux,u_y,a_uy,u_z,a_uz,u_t,a_ut,
v_0,v_x,a_vx,v_y,a_vy,v_z,a_vz,v_t,a_vt,w_0,w_x,a_wx,w_y,a_wy,w_z,a_wz,w_t,a_wt,p_0,p_x,a_px,p_y,a_py,p_z,a_pz,p_t,a_pt,
phi_0,phi_x,a_phix,phi_y,a_phiy,phi_z,a_phiz,phi_t,a_phit,""",real=True)

rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi * t / Lt);
u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt);
v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt);
w_an = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt);
#p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt);
phi_an = phi_0 + phi_x * cos(a_phix * pi * x / L) + phi_y * cos(a_phiy * pi * y / L) + phi_z * sin(a_phiz * pi * z / L) + phi_t * cos(a_phit * pi *
t / Lt);

# Convection: Aplying L1 on rho in order to obtain source term Q1 ----------------------
Q1=L1.subs(u,u_an)
Q1=Q1.subs(v,v_an)
Q1=Q1.subs(w,w_an)
Q1=Q1.subs(rho,rho_an)
Q1=Q1.subs(phi,phi_an)


# Time: Aplying L2 on rho, u, v and w in order to obtain source term Q2 ----------
Q2=L2.subs(u,u_an)
Q2=Q2.subs(v,v_an)
Q2=Q2.subs(w,w_an)
Q2=Q2.subs(rho,rho_an)
Q2=Q2.subs(phi,phi_an)



# Diffusion: Aplying L3 on rho, u, v and w in order to obtain source term Q3 ----------
Q3=L3.subs(u,u_an)
Q3=Q3.subs(v,v_an)
Q3=Q3.subs(w,w_an)
Q3=Q3.subs(rho,rho_an)
Q3=Q3.subs(phi,phi_an)


#--------------------------------------------------------------------------
# Source tem for the continuity equation
Qphi=Q1+Q2+Q3



