
# Loading the file for the equations, manufactured solutions and original soource terms

execfile("navier_stokes_incompressible_g.py")

# Hierarchic MMS ------------------------------------------------------
var('Y,U, V, W, U1,V1,W1')


# Q3 ------------------------------------------------------------------

Q3n=Q3.subs(y**2,Y+1)
Q3n=expand(Q3n)
Q3n=collect(Q3n,[Y*pi**3/(Re*L**3),pi/(Re*L)])


Q3new=Q3n.subs(Y,y**2-1);
R3=Q3-Q3new
R3=expand(R3)
R3==0


# Q2 ------------------------------------------------------------------


L2n=L2.subs(diff(u,x),dudx)
L2n=L2n.subs(diff(u,y),dudy)
L2n=L2n.subs(diff(u,z),dudz)
L2n=L2n.subs(diff(v,x),dvdx)
L2n=L2n.subs(diff(v,y),dvdy)
L2n=L2n.subs(diff(v,z),dvdz)
L2n=L2n.subs(diff(w,x),dwdx)
L2n=L2n.subs(diff(w,y),dwdy)
L2n=L2n.subs(diff(w,z),dwdz)


var('Q2n,U0,Ux,Uy,Uz,Uxy,Uxz,Uyz,V0,Vx,Vy,Vz,Vxy,Vxz,Vyz,W0,Wx,Wy,Wz,Wxy,Wxz,Wyz');


Q2n=Q2.subs(a_u0*cos(g_u0 + f_u0*t) + a_ux*cos(c_ux + b_ux*twopi_invLx*x)*cos(g_ux + f_ux*t) + a_uy*cos(g_uy + f_uy*t)*cos(c_uy + b_uy*twopi_invLy*y) + a_uz*cos(g_uz + f_uz*t)*cos(c_uz + b_uz*twopi_invLz*z) + a_uxy*cos(c_uxy + b_uxy*twopi_invLx*x)*cos(e_uxy + d_uxy*twopi_invLy*y)*cos(g_uxy + f_uxy*t) + a_uxz*cos(g_uxz + f_uxz*t)*cos(c_uxz + b_uxz*twopi_invLx*x)*cos(e_uxz + d_uxz*twopi_invLz*z) + a_uyz*cos(c_uyz + b_uyz*twopi_invLy*y)*cos(g_uyz + f_uyz*t)*cos(e_uyz + d_uyz*twopi_invLz*z),U);
Q2n=Q2.subs(a_v0*cos(g_v0 + f_v0*t) + a_vx*cos(c_vx + b_vx*twopi_invLx*x)*cos(g_vx + f_vx*t) + a_vy*cos(c_vy + b_vy*twopi_invLy*y)*cos(g_vy + f_vy*t) + a_vz*cos(g_vz + f_vz*t)*cos(c_vz + b_vz*twopi_invLz*z) + a_vxy*cos(e_vxy + d_vxy*twopi_invLy*y)*cos(g_vxy + f_vxy*t)*cos(c_vxy + b_vxy*twopi_invLx*x) + a_vxz*cos(e_vxz + d_vxz*twopi_invLz*z)*cos(g_vxz + f_vxz*t)*cos(c_vxz + b_vxz*twopi_invLx*x) + a_vyz*cos(g_vyz + f_vyz*t)*cos(c_vyz + b_vyz*twopi_invLy*y)*cos(e_vyz + d_vyz*twopi_invLz*z),V);

Q2n=Q2.subs(a_w0*cos(g_w0 + f_w0*t) + a_wx*cos(g_wx + f_wx*t)*cos(c_wx + b_wx*twopi_invLx*x) + a_wy*cos(c_wy + b_wy*twopi_invLy*y)*cos(g_wy + f_wy*t) + a_wz*cos(g_wz + f_wz*t)*cos(c_wz + b_wz*twopi_invLz*z) + a_wxy*cos(g_wxy + f_wxy*t)*cos(c_wxy + b_wxy*twopi_invLx*x)*cos(e_wxy + d_wxy*twopi_invLy*y) + a_wxz*cos(g_wxz + f_wxz*t)*cos(c_wxz + b_wxz*twopi_invLx*x)*cos(e_wxz + d_wxz*twopi_invLz*z) + a_wyz*cos(g_wyz + f_wyz*t)*cos(c_wyz + b_wyz*twopi_invLy*y)*cos(e_wyz + d_wyz*twopi_invLz*z),W);  

#Q2n=expand(Q2n)
Q2n=collect(Q2n,[U,V,W])
Q2n=collect([twopi_invLx,twopi_invLy,twopi_invLz])

#Q2n=Q2.subs(  a_u0*cos(g_u0 + f_u0*t),U0);
#Q2n=Q2n.subs( a_ux*cos(c_ux + b_ux*twopi_invLx*x)*cos(g_ux + f_ux*t),Ux);
#Q2n=Q2n.subs( a_uy*cos(g_uy + f_uy*t)*cos(c_uy + b_uy*twopi_invLy*y),Uy); 
#Q2n=Q2n.subs( a_uz*cos(g_uz + f_uz*t)*cos(c_uz + b_uz*twopi_invLz*z),Uz);
#Q2n=Q2n.subs( a_uxy*cos(c_uxy + b_uxy*twopi_invLx*x)*cos(e_uxy + d_uxy*twopi_invLy*y)*cos(g_uxy + f_uxy*t),Uxy);
#Q2n=Q2n.subs( a_uxz*cos(g_uxz + f_uxz*t)*cos(c_uxz + b_uxz*twopi_invLx*x)*cos(e_uxz + d_uxz*twopi_invLz*z),Uxz);
#Q2n=Q2n.subs( a_uyz*cos(c_uyz + b_uyz*twopi_invLy*y)*cos(g_uyz + f_uyz*t)*cos(e_uyz + d_uyz*twopi_invLz*z),Uyz);

#Q2n=Q2n.subs( a_v0*cos(g_v0 + f_v0*t),V0);
#Q2n=Q2n.subs( a_vx*cos(c_vx + b_vx*twopi_invLx*x)*cos(g_vx + f_vx*t) ,Vx);
#Q2n=Q2n.subs( a_vy*cos(c_vy + b_vy*twopi_invLy*y)*cos(g_vy + f_vy*t),Vy);
#Q2n=Q2n.subs( a_vz*cos(g_vz + f_vz*t)*cos(c_vz + b_vz*twopi_invLz*z),Vz);
#Q2n=Q2n.subs( a_vxy*cos(e_vxy + d_vxy*twopi_invLy*y)*cos(g_vxy + f_vxy*t)*cos(c_vxy + b_vxy*twopi_invLx*x),Vxy);
#Q2n=Q2n.subs( a_vxz*cos(e_vxz + d_vxz*twopi_invLz*z)*cos(g_vxz + f_vxz*t)*cos(c_vxz + b_vxz*twopi_invLx*x),Vxz);
#Q2n=Q2n.subs( a_vyz*cos(g_vyz + f_vyz*t)*cos(c_vyz + b_vyz*twopi_invLy*y)*cos(e_vyz + d_vyz*twopi_invLz*z),Vyz);

#Q2n=Q2n.subs(a_w0*cos(g_w0 + f_w0*t) ,W0);
#Q2n=Q2n.subs(a_wx*cos(g_wx + f_wx*t)*cos(c_wx + b_wx*twopi_invLx*x),Wx);
#Q2n=Q2n.subs(a_wy*cos(c_wy + b_wy*twopi_invLy*y)*cos(g_wy + f_wy*t),Wy);
#Q2n=Q2n.subs(a_wz*cos(g_wz + f_wz*t)*cos(c_wz + b_wz*twopi_invLz*z),Wz);
#Q2n=Q2n.subs(a_wxy*cos(g_wxy + f_wxy*t)*cos(c_wxy + b_wxy*twopi_invLx*x)*cos(e_wxy + d_wxy*twopi_invLy*y) ,Wxy);
#Q2n=Q2n.subs(a_wxz*cos(g_wxz + f_wxz*t)*cos(c_wxz + b_wxz*twopi_invLx*x)*cos(e_wxz + d_wxz*twopi_invLz*z) ,Wxz);
#Q2n=Q2n.subs(a_wyz*cos(g_wyz + f_wyz*t)*cos(c_wyz + b_wyz*twopi_invLy*y)*cos(e_wyz + d_wyz*twopi_invLz*z) ,Wyz);

#Q2n=Q2n.subs(U0+Ux+Uy+Uz+Uxy+Uxz+Uyz,U);
#Q2n=Q2n.subs(V0+Vx+Vy+Vz+Vxy+Vxz+Vyz,V);
#Q2n=Q2n.subs(W0+Wx+Wy+Wz+Wxy+Wxz+Wyz,W);

Q2n=expand(Q2n)
Q2n=collect(Q2n,[U,V,W])
Q2n=collect([twopi_invLx,twopi_invLy,twopi_invLz])

Q2n=Q2n.subs(U1*Y,U);
Q2n=Q2n.subs(V1*Y,V);
Q2n=Q2n.subs(W1*Y,W);
Q2n

Q2n=collect(Q2n,[U,V,W]);

Q2n=collect(Q2n,[Y**2*pi**2/L**2])
#pprint(Q2n)



#Checking
Q2new=Q2n.subs(U,(y**2-1)*U1);
Q2new=Q2new.subs(V,(y**2-1)*V1);
Q2new=Q2new.subs(W,(y**2-1)*W1);
Q2new=Q2new.subs(Y,(y**2-1));
Q2new=Q2new.subs(U1,(u_0+u_x*sin(a_ux*pi*x/L)+u_y*cos(a_uy*pi*y/L)+u_z*cos(a_uz*pi*z/L)+u_t*cos(a_ut*pi*t/L)))
Q2new=Q2new.subs(V1,(v_0+v_x*cos(a_vx*pi*x/L)+v_y*sin(a_vy*pi*y/L)+v_z*sin(a_vz*pi*z/L)+v_t*sin(a_vt*pi*t/L)))
Q2new=Q2new.subs(W1,(w_0+w_x*sin(a_wx*pi*x/L)+w_y*sin(a_wy*pi*y/L)+w_z*cos(a_wz*pi*z/L)+w_t*cos(a_wt*pi*t/L)))


Q2new==Q2 # false

Q2new=expand(Q2new)
Q2=expand(Q2)

Q2new==Q2 # true

R2=Q2-Q2new
R2=expand(R2)
R2==0

#---------------------------------------------------
Q_g=Q1+Q2n+Q3n

U=(y**2-1)*(u_0+u_x*sin(a_ux*pi*x/L)+u_y*cos(a_uy*pi*y/L)+u_z*cos(a_uz*pi*z/L)+u_t*cos(a_ut*pi*t/L))
V=(y**2-1)*(v_0+v_x*cos(a_vx*pi*x/L)+v_y*sin(a_vy*pi*y/L)+v_z*sin(a_vz*pi*z/L)+v_t*sin(a_vt*pi*t/L))
W=(y**2-1)*(w_0+w_x*sin(a_wx*pi*x/L)+w_y*sin(a_wy*pi*y/L)+w_z*cos(a_wz*pi*z/L)+w_t*cos(a_wt*pi*t/L))
Y=y**2-1


from sympy.utilities.codegen import codegen,Routine
codegen((
	
	("Y", Y),
	("U", U),
	("V", V),
	("W", W),
        ("Q_g",    Q_g   ),

        ), "C", "Source_Qg", header=True, to_files=True)

