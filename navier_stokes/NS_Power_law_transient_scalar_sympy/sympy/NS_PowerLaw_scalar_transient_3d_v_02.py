
# Loading the file for the equations, manufactured solutions and original soource terms

execfile("NS_PowerLaw_scalar_transient_3d_v_01.py")

# Hierarchic MMS ------------------------------------------------------
var('U,V,W,P,Rho,T')


# Q1 ------------------------------------------------------------------
Q1n=Q1.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi *
t / Lt),Rho);
Q1n=Q1n.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt),U);
Q1n=Q1n.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt),V);
Q1n=Q1n.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt),W);
Q1n=Q1n.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt),P);
Q1n=collect(Q1n, pi*Rho*V/L)


# Q2 ------------------------------------------------------------------
Q2n=Q2



# Q3 ------------------------------------------------------------------
var("Mu, DMu_Dx, DMu_Dy, DMu_Dz,LAMBDA,AUX01")

Q3n=Q3.subs(diff(mu,x),DMu_Dx)
Q3n=Q3n.subs(diff(mu,y),DMu_Dy)
Q3n=Q3n.subs(diff(mu,z),DMu_Dz)
Q3n=Q3n.subs(mu,Mu)
Q3n=Q3n.subs((2/3 - alpha),-LAMBDA/Mu)
Q3n=Q3n.subs(LAMBDA/Mu,AUX01)


Q3n= expand(Q3n) 

Q3n=collect(Q3n,AUX01*DMu_Dy*pi/L)
#Q3n=collect(Q3n,[DMu_Dy*pi/L])
Q3n=collect(Q3n,[DMu_Dz*pi/L])
Q3n=collect(Q3n,[DMu_Dx*pi/L])
Q3n=collect(Q3n,[Mu*pi/L])
Q3n=Q3n.subs(AUX01,alpha-Rational(2,3))


# Q4 ------------------------------------------------------------------
Q4n=Q4.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi *
t / Lt),Rho);
Q4n=Q4n.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt),U);
Q4n=Q4n.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt),V);
Q4n=Q4n.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt),W);
Q4n=Q4n.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt),P);



#--------------------------------------------------------------------------
# Factorized source tem for the u-velocity equation
Q_v_new=Q1n+Q2n+Q3n+Q4n

Q_v_convection=Q1n
Q_v_gradp =Q2n
Q_v_viscous=Q3n
Q_v_time = Q4n


Res=Q_v_new-(Q_v_convection+Q_v_gradp+Q_v_viscous+Q_v_time)
Res==0 #true






##Auxliary calculations for the derivatives of the viscosity mu ->