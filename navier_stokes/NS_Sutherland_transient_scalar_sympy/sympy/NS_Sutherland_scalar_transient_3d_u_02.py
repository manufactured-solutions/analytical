
# Loading the file for the equations, manufactured solutions and original soource terms

execfile("NS_Sutherland_scalar_transient_3d_u_01.py")

# Hierarchic MMS ------------------------------------------------------
var('U,V,W,P,Rho,T')


# Q1 ------------------------------------------------------------------
Q1n=Q1.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi *
t / Lt),Rho);
Q1n=Q1n.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt),U);
Q1n=Q1n.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt),V);
Q1n=Q1n.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt),W);
Q1n=Q1n.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt),P);
Q1n=collect(Q1n, pi*Rho*U/L)


# Q2 ------------------------------------------------------------------
Q2n=Q2


# Q3 ------------------------------------------------------------------
var("Mu, DMu_Dx, DMu_Dy, DMu_Dz")

Q3n=Q3.subs(diff(mu,x),DMu_Dx)
Q3n=Q3n.subs(diff(mu,y),DMu_Dy)
Q3n=Q3n.subs(diff(mu,z),DMu_Dz)
Q3n=Q3n.subs(mu,Mu)

Q3n= expand(Q3n) 

Q3n=collect(Q3n,[DMu_Dy*pi/L])
Q3n=collect(Q3n,[DMu_Dz*pi/L])
Q3n=collect(Q3n,[DMu_Dx*pi/L])
Q3n=collect(Q3n,[Mu*pi/L])

#Q3n=collect(Q3n,[DMu_Dx*pi*L,DMu_Dy*L,DMu_Dz*L,Mu*L])



# Q4 ------------------------------------------------------------------
Q4n=Q4.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi *
t / Lt),Rho);
Q4n=Q4n.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt),U);
Q4n=Q4n.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt),V);
Q4n=Q4n.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt),W);
Q4n=Q4n.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt),P);



#--------------------------------------------------------------------------
# Factorized source tem for the u-velocity equation
Q_u_new=Q1n+Q2n+Q3n+Q4n

Q_u_convection=Q1n
Q_u_gradp =Q2n
Q_u_viscous=Q3n
Q_u_time = Q4n


Res=Q_u_new-(Q_u_convection+Q_u_gradp+Q_u_viscous+Q_u_time)
Res==0 #true






##Auxliary calculations for the derivatives of the viscosity mu
#T = p/(rho*R);
#mu = A_mu*T**(3/2)/(T+B_mu);


## diff(mu,x)
#Dmu_Dx_old=diff(mu,x)
#Dmu_Dx_old=Dmu_Dx_old.subs(u,u_an)
#Dmu_Dx_old=Dmu_Dx_old.subs(v,v_an)
#Dmu_Dx_old=Dmu_Dx_old.subs(w,w_an)
#Dmu_Dx_old=Dmu_Dx_old.subs(rho,rho_an)
#Dmu_Dx_old=Dmu_Dx_old.subs(p,p_an)


## diff(mu,y)
#Dmu_Dy_old=diff(mu,y)
#Dmu_Dy_old=Dmu_Dy_old.subs(u,u_an)
#Dmu_Dy_old=Dmu_Dy_old.subs(v,v_an)
#Dmu_Dy_old=Dmu_Dy_old.subs(w,w_an)
#Dmu_Dy_old=Dmu_Dy_old.subs(rho,rho_an)
#Dmu_Dy_old=Dmu_Dy_old.subs(p,p_an)

## diff(mu,y)
#Dmu_D_oldz=diff(mu,z)
#Dmu_Dz_old=Dmu_Dz_old.subs(u,u_an)
#Dmu_Dz_old=Dmu_Dz_old.subs(v,v_an)
#Dmu_Dz_old=Dmu_Dz_old.subs(w,w_an)
#Dmu_Dz_old=Dmu_Dz_old.subs(rho,rho_an)
#Dmu_Dz_old=Dmu_Dz_old.subs(p,p_an)

## Dmu_Dx ------------------------------------------------------------------
#Dmu_Dx=Dmu_Dx_old.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t *
#sin(a_rhot * pi *
#t / Lt), Rho);
#Dmu_Dx=Dmu_Dx.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t /
#Lt),U);
#Dmu_Dx=Dmu_Dx.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t /
#Lt),V);
#Dmu_Dx=Dmu_Dx.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t /
#Lt),W);
#Dmu_Dx=Dmu_Dx.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t /
#Lt),P);

#Dmu_Dx=Dmu_Dx.subs(P,Rho*R*T)
#Dmu_Dx=sympify(Dmu_Dx)

#print(expand(Dmu_Dx))


## Dmu_Dy ------------------------------------------------------------------
#Dmu_Dy=Dmu_Dy_old.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t *
#sin(a_rhot * pi *
#t / Lt), Rho);
#Dmu_Dy=Dmu_Dy.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t /
#Lt),U);
#Dmu_Dy=Dmu_Dy.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t /
#Lt),V);
#Dmu_Dy=Dmu_Dy.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t /
#Lt),W);
#Dmu_Dy=Dmu_Dy.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t /
#Lt),P);

#Dmu_Dy=Dmu_Dy.subs(P,Rho*R*T)
#Dmu_Dy=sympify(Dmu_Dy)

#print(Dmu_Dy)


## Dmu_Dz ------------------------------------------------------------------
#Dmu_Dz=Dmu_Dz_old.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t *
#sin(a_rhot * pi *
#t / Lt), Rho);
#Dmu_Dz=Dmu_Dz.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t /
#Lt),U);
#Dmu_Dz=Dmu_Dz.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t /
#Lt),V);
#Dmu_Dz=Dmu_Dz.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t /
#Lt),W);
#Dmu_Dz=Dmu_Dz.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t /
#Lt),P);

#Dmu_Dz=Dmu_Dz.subs(P,Rho*R*T)
#Dmu_Dz=sympify(Dmu_Dz)

#print(Dmu_Dz)

