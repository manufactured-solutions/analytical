

#Auxliary calculations for the derivatives of the viscosity mu
T = p/(rho*R);

mu = A_mu*T**(3/2)/(T+B_mu);


# diff(mu,x)
Dmu_Dx_old=diff(mu,x)
Dmu_Dx_old=Dmu_Dx_old.subs(u,u_an)
Dmu_Dx_old=Dmu_Dx_old.subs(v,v_an)
Dmu_Dx_old=Dmu_Dx_old.subs(w,w_an)
Dmu_Dx_old=Dmu_Dx_old.subs(rho,rho_an)
Dmu_Dx_old=Dmu_Dx_old.subs(p,p_an)


# diff(mu,y)
Dmu_Dy_old=diff(mu,y)
Dmu_Dy_old=Dmu_Dy_old.subs(u,u_an)
Dmu_Dy_old=Dmu_Dy_old.subs(v,v_an)
Dmu_Dy_old=Dmu_Dy_old.subs(w,w_an)
Dmu_Dy_old=Dmu_Dy_old.subs(rho,rho_an)
Dmu_Dy_old=Dmu_Dy_old.subs(p,p_an)

# diff(mu,y)
Dmu_Dz_old=diff(mu,z)
Dmu_Dz_old=Dmu_Dz_old.subs(u,u_an)
Dmu_Dz_old=Dmu_Dz_old.subs(v,v_an)
Dmu_Dz_old=Dmu_Dz_old.subs(w,w_an)
Dmu_Dz_old=Dmu_Dz_old.subs(rho,rho_an)
Dmu_Dz_old=Dmu_Dz_old.subs(p,p_an)

var('U,V,W,P,Rho,T')
# Dmu_Dx ------------------------------------------------------------------
Dmu_Dx=Dmu_Dx_old.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t *
sin(a_rhot * pi * t / Lt), Rho);
Dmu_Dx=Dmu_Dx.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t /
Lt),U);
Dmu_Dx=Dmu_Dx.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t /
Lt),V);
Dmu_Dx=Dmu_Dx.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t /
Lt),W);
Dmu_Dx=Dmu_Dx.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t /
Lt),P);

Dmu_Dx=Dmu_Dx.subs(P,Rho*R*T)

Dmu_Dx=sympify(Dmu_Dx)




# Dmu_Dy ------------------------------------------------------------------
Dmu_Dy=Dmu_Dy_old.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t *
sin(a_rhot * pi * t / Lt), Rho);
Dmu_Dy=Dmu_Dy.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t /
Lt),U);
Dmu_Dy=Dmu_Dy.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t /
Lt),V);
Dmu_Dy=Dmu_Dy.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t /
Lt),W);
Dmu_Dy=Dmu_Dy.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t /
Lt),P);

Dmu_Dy=Dmu_Dy.subs(P,Rho*R*T)
Dmu_Dy=sympify(Dmu_Dy)




# Dmu_Dz ------------------------------------------------------------------
Dmu_Dz=Dmu_Dz_old.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t *
sin(a_rhot * pi *
t / Lt), Rho);
Dmu_Dz=Dmu_Dz.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t /
Lt),U);
Dmu_Dz=Dmu_Dz.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t /
Lt),V);
Dmu_Dz=Dmu_Dz.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t /
Lt),W);
Dmu_Dz=Dmu_Dz.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t /
Lt),P);

Dmu_Dz=Dmu_Dz.subs(P,Rho*R*T)
Dmu_Dz=sympify(Dmu_Dz)



# Writing Q_g into a latex file ------------------------------------------



latexQ1=latex(Dmu_Dx)
latexQ2=latex(Dmu_Dy)
latexQ3=latex(Dmu_Dz)

s1=str(latexQ1)
s2=str(latexQ2)
s3=str(latexQ3)

f=open('../latex/derivatives.tex','w')



f.write('\n\n Dmu_Dx: \n')
f.write(s1)

f.write('\n\n Dmu_Dy: \n')
f.write(s2)

f.write('\n\n Dmu_Dz: \n')
f.write(s3)

f.close()


