
# Loading the file for the equations, manufactured solutions and original soource terms

execfile("NS_Sutherland_scalar_transient_3d_phi_01.py")

# Hierarchic MMS ------------------------------------------------------
var('U,V,W,P,Rho,Phi')


# Q1 ------------------------------------------------------------------
Q1n=Q1.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t
* sin(a_rhot * pi * t / Lt),Rho);
Q1n=Q1n.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut *
pi * t / Lt),U);
Q1n=Q1n.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt *
pi * t / Lt),V);
Q1n=Q1n.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt *
pi * t / Lt),W);
Q1n=Q1n.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt *
pi * t / Lt),P);
Q1n=Q1n.subs(phi_0 + phi_x * cos(a_phix * pi * x / L) + phi_y * cos(a_phiy * pi * y / L) + phi_z * sin(a_phiz * pi * z / L) + phi_t * cos(a_phit * pi
*t / Lt) ,Phi);


Q1n=collect(Q1n,[Rho*Phi*pi/L],exact=True)
#Q1n=collect(Q1n,Phi*pi/L)

Q_phi_convection=Q1n

# Q2 ------------------------------------------------------------------
Q2n=Q2.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t
* sin(a_rhot * pi * t / Lt),Rho);
Q2n=Q2n.subs(phi_0 + phi_x * cos(a_phix * pi * x / L) + phi_y * cos(a_phiy * pi * y / L) + phi_z * sin(a_phiz * pi * z / L) + phi_t * cos(a_phit * pi
*t / Lt) ,Phi);

Q_phi_time=Q2n


# Q3 ------------------------------------------------------------------
Q3n=expand(Q3)
Q3n=collect(Q3n,Gamma*pi*L)

Q_phi_diffusion=Q3n



#--------------------------------------------------------------------------
# Checking factorization

Q_phi_new = Q_phi_diffusion +Q_phi_convection + Q_phi_time

Q=Q_phi_new.subs(Rho,rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t *
sin(a_rhot * pi * t / Lt));
Q=Q.subs(U, u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt));
Q=Q.subs(V, v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt));
Q=Q.subs(W, w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt));
Q=Q.subs(P, p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt));
Q=Q.subs(P, phi_0 + phi_x * cos(a_phix * pi * x / L) + phi_y * cos(a_phiy * pi * y / L) + phi_z * sin(a_phiz * pi * z / L) + phi_t * cos(a_phit * pi
*t / Lt));

Res1=expand(Qphi-Q)
Res1==0 #true

# The expression is too complex and Sympy gets lost
Res1=(Qphi-Q);
Res1=Res1.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot *
pi * t / Lt),Rho);
Res1=Res1.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt),U);
Res1=Res1.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt),V);
Res1=Res1.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt),W);
Res1=Res1.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt),P);
Res1=Res1.subs(phi_0 + phi_x * cos(a_phix * pi * x / L) + phi_y * cos(a_phiy * pi * y / L) + phi_z * sin(a_phiz * pi * z / L) + phi_t * cos(a_phit *
pi *t / Lt) ,Phi);
Res1=expand(Res1)
Res1==0 #true



#--------------------------------------------------------------------------
# Writing to C codes
if Res1 == 0:
  print 'Successfull factorization: Q_phi'
  print 'Writing C code for Q_phi -> ../C_codes/NS_Sutherland_scalar_transient_3d_phi.c'
  execfile("NS_Sutherland_scalar_transient_3d_phi_codes.py")  
  print 'Done.'
else:
  print 'ERROR: Possible problems in the factorization.'
  
