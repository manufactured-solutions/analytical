
# Loading the file for the equations, manufactured solutions and original soource terms

execfile("NS_Sutherland_scalar_transient_3d_e_01.py")

# Hierarchic MMS ------------------------------------------------------
var('U,V,W,P,Rho')

var('E_t,G') #Et=e_t

# Q1 ------------------------------------------------------------------
Q1n=Q1.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi *
t / Lt),Rho);
Q1n=Q1n.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt),U);
Q1n=Q1n.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt),V);
Q1n=Q1n.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt),W);
Q1n=Q1n.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt),P);
Q1n=Q1n.subs(-P/(Rho*(1 - gamma)) + U**2/2 + V**2/2 + W**2/2,E_t)
#Q1n=Q1n.subs(1 - gamma, G)
A=expand(Q1n)
A=collect(A,[E_t*Rho*L*pi])

A1=pi*E_t*Rho*(a_ux*u_x*cos(pi*a_ux*x/L) + a_vy*v_y*cos(pi*a_vy*y/L) - a_wz*w_z*sin(pi*a_wz*z/L))/L

A2=A-A1
A2= collect(A2,[pi*a_rhox*rho_x*cos(pi*a_rhox*x/L)*L])
A2= collect(A2,[pi*a_rhoz*rho_z*cos(pi*a_rhoz*z/L)*L])
A2= collect(A2,[pi*a_rhoy*rho_y*sin(pi*a_rhoy*y/L)*L])


A3=pi*a_rhox*rho_x*(E_t*U + P*U/(Rho*(1 - gamma)))*cos(pi*a_rhox*x/L)/L + pi*a_rhoy*rho_y*(-E_t*V - P*V/(Rho*(1
-gamma)))*sin(pi*a_rhoy*y/L)/L+pi*a_rhoz*rho_z*(E_t*W + P*W/(Rho*(1 - gamma)))*cos(pi*a_rhoz*z/L)/L 

A4=A2-A3
A4=collect(A4,[Rho*L*pi])

A4=collect(A4,(1-gamma)*L*pi)


A2=A3+A4
A=A1+A2
Q1n=A


###print(Q1n)
##A1=Rho*U*(pi*U*a_ux*u_x*cos(pi*a_ux*x/L)/L + pi*W*a_wx*w_x*cos(pi*a_wx*x/L)/L
#-pi*V*a_vx*v_x*sin(pi*a_vx*x/L)/L+pi*a_px*p_x*sin(pi*a_px*x/L)/(L*Rho*(1 - gamma)) + pi*P*a_rhox*rho_x*cos(pi*a_rhox*x/L)/(L*Rho**2*(1 - gamma)))
#+Rho*V*(pi*V*a_vy*v_y*cos(pi*a_vy*y/L)/L + pi*W*a_wy*w_y*cos(pi*a_wy*y/L)/L - pi*U*a_uy*u_y*sin(pi*a_uy*y/L)/L
#-pi*a_py*p_y*cos(pi*a_py*y/L)/(L*Rho*(1- gamma)) - pi*P*a_rhoy*rho_y*sin(pi*a_rhoy*y/L)/(L*Rho**2*(1 - gamma)))
#+Rho*W*(pi*V*a_vz*v_z*cos(pi*a_vz*z/L)/L - pi*U*a_uz*u_z*sin(pi*a_uz*z/L)/L - pi*W*a_wz*w_z*sin(pi*a_wz*z/L)/L
#+pi*a_pz*p_z*sin(pi*a_pz*z/L)/(L*Rho*(1- gamma)) + pi*P*a_rhoz*rho_z*cos(pi*a_rhoz*z/L)/(L*Rho**2*(1 - gamma))) 

#A2=Q1n-A1;
#A2=collect(A2,[pi*(-P/(Rho*(1 - gamma)) + U**2/2 + V**2/2 + W**2/2)/L])

#Q1nn=A1+A2;


# Q2 ------------------------------------------------------------------
Q2n=Q2.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi *
t / Lt),Rho);
Q2n=Q2n.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt),U);
Q2n=Q2n.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt),V);
Q2n=Q2n.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt),W);
Q2n=Q2n.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt),P);

Q2n=collect(Q2n,P*pi/L)


# Q3 ------------------------------------------------------------------
var("Mu, DMu_Dx, DMu_Dy, DMu_Dz")
Q3n=Q3.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi *
t / Lt),Rho);
Q3n=Q3n.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt),U);
Q3n=Q3n.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt),V);
Q3n=Q3n.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt),W);
Q3n=Q3n.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt),P);
Q3n=Q3n.subs(diff(mu,x),DMu_Dx)
Q3n=Q3n.subs(diff(mu,y),DMu_Dy)
Q3n=Q3n.subs(diff(mu,z),DMu_Dz)
Q3n=Q3n.subs(mu,Mu)
Q3n= expand(Q3n) 


Q3n=collect(Q3n,[U*DMu_Dx*pi/L])
Q3n=collect(Q3n,[V*DMu_Dx*pi/L])
Q3n=collect(Q3n,[W*DMu_Dx*pi/L])

Q3n=collect(Q3n,[U*DMu_Dy*pi/L])
Q3n=collect(Q3n,[V*DMu_Dy*pi/L])
Q3n=collect(Q3n,[W*DMu_Dy*pi/L])

Q3n=collect(Q3n,[U*DMu_Dz*pi/L])
Q3n=collect(Q3n,[V*DMu_Dz*pi/L])
Q3n=collect(Q3n,[W*DMu_Dz*pi/L])


Q3n=collect(Q3n,[U*Mu*pi/L])
Q3n=collect(Q3n,[V*Mu*pi/L])
Q3n=collect(Q3n,[W*Mu*pi/L])

Q3n=collect(Q3n,[Mu*pi/L])



# Q4 ------------------------------------------------------------------
Q4n=Q4.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi *
t / Lt),Rho);
Q4n=Q4n.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt),U);
Q4n=Q4n.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt),V);
Q4n=Q4n.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt),W);
Q4n=Q4n.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt),P);

Q4n=expand(Q4n)
Q4n=collect(Q4n,[2*k*pi**2*P/(L**2*R*Rho**3)],exact=True)
Q4n=collect(Q4n,[k*pi**2*P/(L**2*R*Rho**2)],exact=True)

#print(Q4n)
A1=P*k*pi**2*(-2*a_rhox**2*rho_x**2*cos(pi*a_rhox*x/L)**2 - 2*a_rhoy**2*rho_y**2*sin(pi*a_rhoy*y/L)**2 - 2*a_rhoz**2*rho_z**2*cos(pi*a_rhoz*z/L)**2)/(L**2*R*Rho**3) + P*k*pi**2*(-rho_x*a_rhox**2*sin(pi*a_rhox*x/L) - rho_y*a_rhoy**2*cos(pi*a_rhoy*y/L) - rho_z*a_rhoz**2*sin(pi*a_rhoz*z/L))/(L**2*R*Rho**2) 
A2=Q4n-A1

A2=collect(A2,[k*pi**2/(L**2*R*Rho)],exact=True)
A2=collect(A2,[k*pi**2/(L**2*R*Rho**2)],exact=True)

Q4new=A1+A2

Res=expand(Q4n-Q4new)
Res==0 #True

Q4n=Q4new


# Q5 ------------------------------------------------------------------
Q5n=Q5.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi *
t / Lt),Rho);
Q5n=Q5n.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt),U);
Q5n=Q5n.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt),V);
Q5n=Q5n.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt),W);
Q5n=Q5n.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt),P);

Q5n=Q5n.subs(-P/(Rho*(1 - gamma)) + U**2/2 + V**2/2 + W**2/2,E_t)
Q5n=expand(Q5n)
Q5n=collect(Q5n,Rho*Lt*pi)




#--------------------------------------------------------------------------
# Factorized source tem for the u-velocity equation
Q_et_new=Q1n+Q2n+Q3n+Q4n+Q5n

Q_et_convection=Q1n
Q_et_gradp =Q2n
Q_et_viscous=Q3n
Q_et_heatflux=Q4n
Q_et_time = Q5n


Res=Q_et_new-(Q_et_convection+Q_et_gradp+Q_et_viscous+Q_et_heatflux+Q_et_time)
Res==0 #true

