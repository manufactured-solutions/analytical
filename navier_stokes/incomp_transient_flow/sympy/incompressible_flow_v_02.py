
# Loading the file for the equations, manufactured solutions and original soource terms

execfile("incompressible_flow_v.py")

# Hierarchic MMS ------------------------------------------------------
var('Y,U, V, W, Q1n,Q2n,Q3n,R1,R2,R3')
# Q1 ------------------------------------------------------------------
Q1n=collect(Q1,[twopi_invLx,twopi_invLz])

R1=Q1-Q1n
R1=expand(R1)
R1==0 #true

# Q3 ------------------------------------------------------------------

Q3n=expand(Q3)

Q3n=collect(Q3n,[twopi_invLx**4/Re],exact=True)
Q3n=collect(Q3n,[twopi_invLy**4/Re],exact=True)
Q3n=collect(Q3n,[twopi_invLz**4/Re],exact=True)

R3=Q3-Q3n
R3=expand(R3)
R3==0 #true



# Q2 ------------------------------------------------------------------
var('U0,Ux,Uy,Uz,Uxy,Uxz,Uyz,V0,Vx,Vy,Vz,Vxy,Vxz,Vyz,W0,Wx,Wy,Wz,Wxy,Wxz,Wyz');

Q2n=Q2.subs(  a_u0*cos(g_u0 + f_u0*t),U0);
Q2n=Q2n.subs( a_ux*cos(c_ux + b_ux*twopi_invLx*x)*cos(g_ux + f_ux*t),Ux);
Q2n=Q2n.subs( a_uy*cos(g_uy + f_uy*t)*cos(c_uy + b_uy*twopi_invLy*y),Uy); 
Q2n=Q2n.subs( a_uz*cos(g_uz + f_uz*t)*cos(c_uz + b_uz*twopi_invLz*z),Uz);
Q2n=Q2n.subs( a_uxy*cos(c_uxy + b_uxy*twopi_invLx*x)*cos(e_uxy + d_uxy*twopi_invLy*y)*cos(g_uxy + f_uxy*t),Uxy);
Q2n=Q2n.subs( a_uxz*cos(g_uxz + f_uxz*t)*cos(c_uxz + b_uxz*twopi_invLx*x)*cos(e_uxz + d_uxz*twopi_invLz*z),Uxz);
Q2n=Q2n.subs( a_uyz*cos(c_uyz + b_uyz*twopi_invLy*y)*cos(g_uyz + f_uyz*t)*cos(e_uyz + d_uyz*twopi_invLz*z),Uyz);

Q2n=Q2n.subs( a_v0*cos(g_v0 + f_v0*t),V0);
Q2n=Q2n.subs( a_vx*cos(c_vx + b_vx*twopi_invLx*x)*cos(g_vx + f_vx*t) ,Vx);
Q2n=Q2n.subs( a_vy*cos(c_vy + b_vy*twopi_invLy*y)*cos(g_vy + f_vy*t),Vy);
Q2n=Q2n.subs( a_vz*cos(g_vz + f_vz*t)*cos(c_vz + b_vz*twopi_invLz*z),Vz);
Q2n=Q2n.subs( a_vxy*cos(e_vxy + d_vxy*twopi_invLy*y)*cos(g_vxy + f_vxy*t)*cos(c_vxy + b_vxy*twopi_invLx*x),Vxy);
Q2n=Q2n.subs( a_vxz*cos(e_vxz + d_vxz*twopi_invLz*z)*cos(g_vxz + f_vxz*t)*cos(c_vxz + b_vxz*twopi_invLx*x),Vxz);
Q2n=Q2n.subs( a_vyz*cos(g_vyz + f_vyz*t)*cos(c_vyz + b_vyz*twopi_invLy*y)*cos(e_vyz + d_vyz*twopi_invLz*z),Vyz);

Q2n=Q2n.subs(a_w0*cos(g_w0 + f_w0*t) ,W0);
Q2n=Q2n.subs(a_wx*cos(g_wx + f_wx*t)*cos(c_wx + b_wx*twopi_invLx*x),Wx);
Q2n=Q2n.subs(a_wy*cos(c_wy + b_wy*twopi_invLy*y)*cos(g_wy + f_wy*t),Wy);
Q2n=Q2n.subs(a_wz*cos(g_wz + f_wz*t)*cos(c_wz + b_wz*twopi_invLz*z),Wz);
Q2n=Q2n.subs(a_wxy*cos(g_wxy + f_wxy*t)*cos(c_wxy + b_wxy*twopi_invLx*x)*cos(e_wxy + d_wxy*twopi_invLy*y) ,Wxy);
Q2n=Q2n.subs(a_wxz*cos(g_wxz + f_wxz*t)*cos(c_wxz + b_wxz*twopi_invLx*x)*cos(e_wxz + d_wxz*twopi_invLz*z) ,Wxz);
Q2n=Q2n.subs(a_wyz*cos(g_wyz + f_wyz*t)*cos(c_wyz + b_wyz*twopi_invLy*y)*cos(e_wyz + d_wyz*twopi_invLz*z) ,Wyz);

Q2n=Q2n.subs(U0+Ux+Uy+Uz+Uxy+Uxz+Uyz,U);
Q2n=Q2n.subs(V0+Vx+Vy+Vz+Vxy+Vxz+Vyz,V);
Q2n=Q2n.subs(W0+Wx+Wy+Wz+Wxy+Wxz+Wyz,W);

Q2n=collect(Q2n,[U,V,W])
#print(Q2n)
 
A1=U*(-a_vx*b_vx**3*twopi_invLx**3*cos(g_vx + f_vx*t)*sin(c_vx + b_vx*twopi_invLx*x) - a_vxy*b_vxy**3*twopi_invLx**3*cos(e_vxy + d_vxy*twopi_invLy*y)*cos(g_vxy +f_vxy*t)*sin(c_vxy + b_vxy*twopi_invLx*x) - a_vxz*b_vxz**3*twopi_invLx**3*cos(e_vxz + d_vxz*twopi_invLz*z)*cos(g_vxz + f_vxz*t)*sin(c_vxz + b_vxz*twopi_invLx*x)+ a_uxy*d_uxy*twopi_invLy*b_uxy**2*twopi_invLx**2*cos(c_uxy + b_uxy*twopi_invLx*x)*cos(g_uxy + f_uxy*t)*sin(e_uxy + d_uxy*twopi_invLy*y) -a_vxz*b_vxz*twopi_invLx*d_vxz**2*twopi_invLz**2*cos(e_vxz + d_vxz*twopi_invLz*z)*cos(g_vxz + f_vxz*t)*sin(c_vxz + b_vxz*twopi_invLx*x)) +V*(a_uxy*b_uxy*twopi_invLx*d_uxy**2*twopi_invLy**2*cos(e_uxy + d_uxy*twopi_invLy*y)*cos(g_uxy + f_uxy*t)*sin(c_uxy + b_uxy*twopi_invLx*x) +a_wyz*d_wyz*twopi_invLz*b_wyz**2*twopi_invLy**2*cos(g_wyz + f_wyz*t)*cos(c_wyz + b_wyz*twopi_invLy*y)*sin(e_wyz + d_wyz*twopi_invLz*z) -a_vxy*d_vxy*twopi_invLy*b_vxy**2*twopi_invLx**2*cos(g_vxy + f_vxy*t)*cos(c_vxy + b_vxy*twopi_invLx*x)*sin(e_vxy + d_vxy*twopi_invLy*y) -a_vyz*b_vyz*twopi_invLy*d_vyz**2*twopi_invLz**2*cos(g_vyz + f_vyz*t)*cos(e_vyz + d_vyz*twopi_invLz*z)*sin(c_vyz + b_vyz*twopi_invLy*y)) +W*(-a_vz*b_vz**3*twopi_invLz**3*cos(g_vz + f_vz*t)*sin(c_vz + b_vz*twopi_invLz*z) - a_vxz*d_vxz**3*twopi_invLz**3*cos(g_vxz + f_vxz*t)*cos(c_vxz +b_vxz*twopi_invLx*x)*sin(e_vxz + d_vxz*twopi_invLz*z) - a_vyz*d_vyz**3*twopi_invLz**3*cos(g_vyz + f_vyz*t)*cos(c_vyz + b_vyz*twopi_invLy*y)*sin(e_vyz +d_vyz*twopi_invLz*z) + a_wyz*b_wyz*twopi_invLy*d_wyz**2*twopi_invLz**2*cos(g_wyz + f_wyz*t)*cos(e_wyz + d_wyz*twopi_invLz*z)*sin(c_wyz + b_wyz*twopi_invLy*y) -a_vxz*d_vxz*twopi_invLz*b_vxz**2*twopi_invLx**2*cos(g_vxz + f_vxz*t)*cos(c_vxz + b_vxz*twopi_invLx*x)*sin(e_vxz + d_vxz*twopi_invLz*z))


A2=Q2n-A1
#A2=(-Ux*b_ux**2*twopi_invLx**2 - Uxy*b_uxy**2*twopi_invLx**2 - Uxz*b_uxz**2*twopi_invLx**2)*(-a_uy*b_uy*twopi_invLy*cos(g_uy + f_uy*t)*sin(c_uy + b_uy*twopi_invLy*y) - a_uxy*d_uxy*twopi_invLy*cos(c_uxy + b_uxy*twopi_invLx*x)*cos(g_uxy + f_uxy*t)*sin(e_uxy + d_uxy*twopi_invLy*y) - a_uyz*b_uyz*twopi_invLy*cos(g_uyz + f_uyz*t)*cos(e_uyz + d_uyz*twopi_invLz*z)*sin(c_uyz + b_uyz*twopi_invLy*y)) + (-Uxy*d_uxy**2*twopi_invLy**2 - Uy*b_uy**2*twopi_invLy**2 - Uyz*b_uyz**2*twopi_invLy**2)*(-a_vx*b_vx*twopi_invLx*cos(g_vx + f_vx*t)*sin(c_vx + b_vx*twopi_invLx*x) - a_vxy*b_vxy*twopi_invLx*cos(e_vxy + d_vxy*twopi_invLy*y)*cos(g_vxy + f_vxy*t)*sin(c_vxy + b_vxy*twopi_invLx*x) - a_vxz*b_vxz*twopi_invLx*cos(e_vxz + d_vxz*twopi_invLz*z)*cos(g_vxz + f_vxz*t)*sin(c_vxz + b_vxz*twopi_invLx*x)) + (-Wxy*d_wxy**2*twopi_invLy**2 - Wy*b_wy**2*twopi_invLy**2 - Wyz*b_wyz**2*twopi_invLy**2)*(-a_vz*b_vz*twopi_invLz*cos(g_vz + f_vz*t)*sin(c_vz + b_vz*twopi_invLz*z) - a_vxz*d_vxz*twopi_invLz*cos(g_vxz + f_vxz*t)*cos(c_vxz + b_vxz*twopi_invLx*x)*sin(e_vxz + d_vxz*twopi_invLz*z) - a_vyz*d_vyz*twopi_invLz*cos(g_vyz + f_vyz*t)*cos(c_vyz + b_vyz*twopi_invLy*y)*sin(e_vyz + d_vyz*twopi_invLz*z)) + (-Wxz*d_wxz**2*twopi_invLz**2 - Wyz*d_wyz**2*twopi_invLz**2 - Wz*b_wz**2*twopi_invLz**2)*(-a_wy*b_wy*twopi_invLy*cos(g_wy + f_wy*t)*sin(c_wy + b_wy*twopi_invLy*y) - a_wxy*d_wxy*twopi_invLy*cos(g_wxy + f_wxy*t)*cos(c_wxy + b_wxy*twopi_invLx*x)*sin(e_wxy + d_wxy*twopi_invLy*y) - a_wyz*b_wyz*twopi_invLy*cos(g_wyz + f_wyz*t)*cos(e_wyz + d_wyz*twopi_invLz*z)*sin(c_wyz + b_wyz*twopi_invLy*y)) - (-Ux*b_ux**2*twopi_invLx**2 - Uxy*b_uxy**2*twopi_invLx**2 - Uxz*b_uxz**2*twopi_invLx**2)*(-a_vx*b_vx*twopi_invLx*cos(g_vx + f_vx*t)*sin(c_vx + b_vx*twopi_invLx*x) - a_vxy*b_vxy*twopi_invLx*cos(e_vxy + d_vxy*twopi_invLy*y)*cos(g_vxy + f_vxy*t)*sin(c_vxy + b_vxy*twopi_invLx*x) - a_vxz*b_vxz*twopi_invLx*cos(e_vxz + d_vxz*twopi_invLz*z)*cos(g_vxz + f_vxz*t)*sin(c_vxz + b_vxz*twopi_invLx*x)) - (-Uxz*d_uxz**2*twopi_invLz**2 - Uyz*d_uyz**2*twopi_invLz**2 - Uz*b_uz**2*twopi_invLz**2)*(-a_vx*b_vx*twopi_invLx*cos(g_vx + f_vx*t)*sin(c_vx + b_vx*twopi_invLx*x) - a_vxy*b_vxy*twopi_invLx*cos(e_vxy + d_vxy*twopi_invLy*y)*cos(g_vxy + f_vxy*t)*sin(c_vxy + b_vxy*twopi_invLx*x) - a_vxz*b_vxz*twopi_invLx*cos(e_vxz + d_vxz*twopi_invLz*z)*cos(g_vxz + f_vxz*t)*sin(c_vxz + b_vxz*twopi_invLx*x)) - (-Vx*b_vx**2*twopi_invLx**2 - Vxy*b_vxy**2*twopi_invLx**2 - Vxz*b_vxz**2*twopi_invLx**2)*(-a_vy*b_vy*twopi_invLy*cos(g_vy + f_vy*t)*sin(c_vy + b_vy*twopi_invLy*y) - a_vxy*d_vxy*twopi_invLy*cos(g_vxy + f_vxy*t)*cos(c_vxy + b_vxy*twopi_invLx*x)*sin(e_vxy + d_vxy*twopi_invLy*y) - a_vyz*b_vyz*twopi_invLy*cos(g_vyz + f_vyz*t)*cos(e_vyz + d_vyz*twopi_invLz*z)*sin(c_vyz + b_vyz*twopi_invLy*y)) - (-Vxz*d_vxz**2*twopi_invLz**2 - Vyz*d_vyz**2*twopi_invLz**2 - Vz*b_vz**2*twopi_invLz**2)*(-a_vy*b_vy*twopi_invLy*cos(g_vy + f_vy*t)*sin(c_vy + b_vy*twopi_invLy*y) - a_vxy*d_vxy*twopi_invLy*cos(g_vxy + f_vxy*t)*cos(c_vxy + b_vxy*twopi_invLx*x)*sin(e_vxy + d_vxy*twopi_invLy*y) - a_vyz*b_vyz*twopi_invLy*cos(g_vyz + f_vyz*t)*cos(e_vyz + d_vyz*twopi_invLz*z)*sin(c_vyz + b_vyz*twopi_invLy*y)) - (-Wx*b_wx**2*twopi_invLx**2 - Wxy*b_wxy**2*twopi_invLx**2 - Wxz*b_wxz**2*twopi_invLx**2)*(-a_vz*b_vz*twopi_invLz*cos(g_vz + f_vz*t)*sin(c_vz + b_vz*twopi_invLz*z) - a_vxz*d_vxz*twopi_invLz*cos(g_vxz + f_vxz*t)*cos(c_vxz + b_vxz*twopi_invLx*x)*sin(e_vxz + d_vxz*twopi_invLz*z) - a_vyz*d_vyz*twopi_invLz*cos(g_vyz + f_vyz*t)*cos(c_vyz + b_vyz*twopi_invLy*y)*sin(e_vyz + d_vyz*twopi_invLz*z)) - (-Wxz*d_wxz**2*twopi_invLz**2 - Wyz*d_wyz**2*twopi_invLz**2 - Wz*b_wz**2*twopi_invLz**2)*(-a_vz*b_vz*twopi_invLz*cos(g_vz + f_vz*t)*sin(c_vz + b_vz*twopi_invLz*z) - a_vxz*d_vxz*twopi_invLz*cos(g_vxz + f_vxz*t)*cos(c_vxz + b_vxz*twopi_invLx*x)*sin(e_vxz + d_vxz*twopi_invLz*z) - a_vyz*d_vyz*twopi_invLz*cos(g_vyz + f_vyz*t)*cos(c_vyz + b_vyz*twopi_invLy*y)*sin(e_vyz + d_vyz*twopi_invLz*z)) - 2*(-Vx*b_vx**2*twopi_invLx**2 - Vxy*b_vxy**2*twopi_invLx**2 - Vxz*b_vxz**2*twopi_invLx**2)*(-a_ux*b_ux*twopi_invLx*cos(g_ux + f_ux*t)*sin(c_ux + b_ux*twopi_invLx*x) - a_uxy*b_uxy*twopi_invLx*cos(e_uxy + d_uxy*twopi_invLy*y)*cos(g_uxy + f_uxy*t)*sin(c_uxy + b_uxy*twopi_invLx*x) - a_uxz*b_uxz*twopi_invLx*cos(g_uxz + f_uxz*t)*cos(e_uxz + d_uxz*twopi_invLz*z)*sin(c_uxz + b_uxz*twopi_invLx*x)) - 2*(-Vxz*d_vxz**2*twopi_invLz**2 - Vyz*d_vyz**2*twopi_invLz**2 - Vz*b_vz**2*twopi_invLz**2)*(-a_wz*b_wz*twopi_invLz*cos(g_wz + f_wz*t)*sin(c_wz + b_wz*twopi_invLz*z) - a_wxz*d_wxz*twopi_invLz*cos(g_wxz + f_wxz*t)*cos(c_wxz + b_wxz*twopi_invLx*x)*sin(e_wxz + d_wxz*twopi_invLz*z) - a_wyz*d_wyz*twopi_invLz*cos(g_wyz + f_wyz*t)*cos(c_wyz + b_wyz*twopi_invLy*y)*sin(e_wyz + d_wyz*twopi_invLz*z))+a_uxy*b_uxy*d_uxy*twopi_invLx*twopi_invLy*(-a_vy*b_vy*twopi_invLy*cos(g_vy + f_vy*t)*sin(c_vy + b_vy*twopi_invLy*y) - a_vxy*d_vxy*twopi_invLy*cos(g_vxy + f_vxy*t)*cos(c_vxy + b_vxy*twopi_invLx*x)*sin(e_vxy + d_vxy*twopi_invLy*y) - a_vyz*b_vyz*twopi_invLy*cos(g_vyz + f_vyz*t)*cos(e_vyz + d_vyz*twopi_invLz*z)*sin(c_vyz + b_vyz*twopi_invLy*y))*cos(g_uxy + f_uxy*t)*sin(c_uxy + b_uxy*twopi_invLx*x)*sin(e_uxy + d_uxy*twopi_invLy*y) + a_uxz*b_uxz*d_uxz*twopi_invLx*twopi_invLz*(-a_wy*b_wy*twopi_invLy*cos(g_wy + f_wy*t)*sin(c_wy + b_wy*twopi_invLy*y) - a_wxy*d_wxy*twopi_invLy*cos(g_wxy + f_wxy*t)*cos(c_wxy + b_wxy*twopi_invLx*x)*sin(e_wxy + d_wxy*twopi_invLy*y) - a_wyz*b_wyz*twopi_invLy*cos(g_wyz + f_wyz*t)*cos(e_wyz + d_wyz*twopi_invLz*z)*sin(c_wyz + b_wyz*twopi_invLy*y))*cos(g_uxz + f_uxz*t)*sin(c_uxz + b_uxz*twopi_invLx*x)*sin(e_uxz + d_uxz*twopi_invLz*z) + a_vxy*b_vxy*d_vxy*twopi_invLx*twopi_invLy*(-a_uy*b_uy*twopi_invLy*cos(g_uy + f_uy*t)*sin(c_uy + b_uy*twopi_invLy*y) - a_uxy*d_uxy*twopi_invLy*cos(c_uxy + b_uxy*twopi_invLx*x)*cos(g_uxy + f_uxy*t)*sin(e_uxy + d_uxy*twopi_invLy*y) - a_uyz*b_uyz*twopi_invLy*cos(g_uyz + f_uyz*t)*cos(e_uyz + d_uyz*twopi_invLz*z)*sin(c_uyz + b_uyz*twopi_invLy*y))*cos(g_vxy + f_vxy*t)*sin(e_vxy + d_vxy*twopi_invLy*y)*sin(c_vxy + b_vxy*twopi_invLx*x) + a_vyz*b_vyz*d_vyz*twopi_invLy*twopi_invLz*(-a_wy*b_wy*twopi_invLy*cos(g_wy + f_wy*t)*sin(c_wy + b_wy*twopi_invLy*y) - a_wxy*d_wxy*twopi_invLy*cos(g_wxy + f_wxy*t)*cos(c_wxy + b_wxy*twopi_invLx*x)*sin(e_wxy + d_wxy*twopi_invLy*y) - a_wyz*b_wyz*twopi_invLy*cos(g_wyz + f_wyz*t)*cos(e_wyz + d_wyz*twopi_invLz*z)*sin(c_wyz + b_wyz*twopi_invLy*y))*cos(g_vyz + f_vyz*t)*sin(c_vyz + b_vyz*twopi_invLy*y)*sin(e_vyz + d_vyz*twopi_invLz*z) + a_wxz*b_wxz*d_wxz*twopi_invLx*twopi_invLz*(-a_uy*b_uy*twopi_invLy*cos(g_uy + f_uy*t)*sin(c_uy + b_uy*twopi_invLy*y) - a_uxy*d_uxy*twopi_invLy*cos(c_uxy + b_uxy*twopi_invLx*x)*cos(g_uxy + f_uxy*t)*sin(e_uxy + d_uxy*twopi_invLy*y) - a_uyz*b_uyz*twopi_invLy*cos(g_uyz + f_uyz*t)*cos(e_uyz + d_uyz*twopi_invLz*z)*sin(c_uyz + b_uyz*twopi_invLy*y))*cos(g_wxz + f_wxz*t)*sin(c_wxz + b_wxz*twopi_invLx*x)*sin(e_wxz + d_wxz*twopi_invLz*z) + a_wyz*b_wyz*d_wyz*twopi_invLy*twopi_invLz*(-a_vy*b_vy*twopi_invLy*cos(g_vy + f_vy*t)*sin(c_vy + b_vy*twopi_invLy*y) - a_vxy*d_vxy*twopi_invLy*cos(g_vxy + f_vxy*t)*cos(c_vxy + b_vxy*twopi_invLx*x)*sin(e_vxy + d_vxy*twopi_invLy*y) - a_vyz*b_vyz*twopi_invLy*cos(g_vyz + f_vyz*t)*cos(e_vyz + d_vyz*twopi_invLz*z)*sin(c_vyz + b_vyz*twopi_invLy*y))*cos(g_wyz + f_wyz*t)*sin(c_wyz + b_wyz*twopi_invLy*y)*sin(e_wyz + d_wyz*twopi_invLz*z) - 2*a_vxy*b_vxy*d_vxy*twopi_invLx*twopi_invLy*(-a_vx*b_vx*twopi_invLx*cos(g_vx + f_vx*t)*sin(c_vx + b_vx*twopi_invLx*x) - a_vxy*b_vxy*twopi_invLx*cos(e_vxy + d_vxy*twopi_invLy*y)*cos(g_vxy + f_vxy*t)*sin(c_vxy + b_vxy*twopi_invLx*x) - a_vxz*b_vxz*twopi_invLx*cos(e_vxz + d_vxz*twopi_invLz*z)*cos(g_vxz + f_vxz*t)*sin(c_vxz + b_vxz*twopi_invLx*x))*cos(g_vxy + f_vxy*t)*sin(e_vxy + d_vxy*twopi_invLy*y)*sin(c_vxy + b_vxy*twopi_invLx*x) - 2*a_vxz*b_vxz*d_vxz*twopi_invLx*twopi_invLz*(-a_wx*b_wx*twopi_invLx*cos(g_wx + f_wx*t)*sin(c_wx + b_wx*twopi_invLx*x) - a_wxy*b_wxy*twopi_invLx*cos(g_wxy + f_wxy*t)*cos(e_wxy + d_wxy*twopi_invLy*y)*sin(c_wxy + b_wxy*twopi_invLx*x) - a_wxz*b_wxz*twopi_invLx*cos(g_wxz + f_wxz*t)*cos(e_wxz + d_wxz*twopi_invLz*z)*sin(c_wxz + b_wxz*twopi_invLx*x))*cos(g_vxz + f_vxz*t)*sin(e_vxz + d_vxz*twopi_invLz*z)*sin(c_vxz + b_vxz*twopi_invLx*x) - 2*a_vxz*b_vxz*d_vxz*twopi_invLx*twopi_invLz*(-a_uz*b_uz*twopi_invLz*cos(g_uz + f_uz*t)*sin(c_uz + b_uz*twopi_invLz*z) - a_uxz*d_uxz*twopi_invLz*cos(g_uxz + f_uxz*t)*cos(c_uxz + b_uxz*twopi_invLx*x)*sin(e_uxz + d_uxz*twopi_invLz*z) - a_uyz*d_uyz*twopi_invLz*cos(c_uyz + b_uyz*twopi_invLy*y)*cos(g_uyz + f_uyz*t)*sin(e_uyz + d_uyz*twopi_invLz*z))*cos(g_vxz + f_vxz*t)*sin(e_vxz + d_vxz*twopi_invLz*z)*sin(c_vxz + b_vxz*twopi_invLx*x) - 2*a_vyz*b_vyz*d_vyz*twopi_invLy*twopi_invLz*(-a_vz*b_vz*twopi_invLz*cos(g_vz + f_vz*t)*sin(c_vz + b_vz*twopi_invLz*z) - a_vxz*d_vxz*twopi_invLz*cos(g_vxz + f_vxz*t)*cos(c_vxz + b_vxz*twopi_invLx*x)*sin(e_vxz + d_vxz*twopi_invLz*z) - a_vyz*d_vyz*twopi_invLz*cos(g_vyz + f_vyz*t)*cos(c_vyz + b_vyz*twopi_invLy*y)*sin(e_vyz + d_vyz*twopi_invLz*z))*cos(g_vyz + f_vyz*t)*sin(c_vyz + b_vyz*twopi_invLy*y)*sin(e_vyz + d_vyz*twopi_invLz*z) + 2*a_uxy*b_uxy*d_uxy*twopi_invLx*twopi_invLy*(-a_ux*b_ux*twopi_invLx*cos(g_ux + f_ux*t)*sin(c_ux + b_ux*twopi_invLx*x) - a_uxy*b_uxy*twopi_invLx*cos(e_uxy + d_uxy*twopi_invLy*y)*cos(g_uxy + f_uxy*t)*sin(c_uxy + b_uxy*twopi_invLx*x) - a_uxz*b_uxz*twopi_invLx*cos(g_uxz + f_uxz*t)*cos(e_uxz + d_uxz*twopi_invLz*z)*sin(c_uxz + b_uxz*twopi_invLx*x))*cos(g_uxy + f_uxy*t)*sin(c_uxy + b_uxy*twopi_invLx*x)*sin(e_uxy + d_uxy*twopi_invLy*y) + 2*a_uyz*b_uyz*d_uyz*twopi_invLy*twopi_invLz*(-a_wx*b_wx*twopi_invLx*cos(g_wx + f_wx*t)*sin(c_wx + b_wx*twopi_invLx*x) - a_wxy*b_wxy*twopi_invLx*cos(g_wxy + f_wxy*t)*cos(e_wxy + d_wxy*twopi_invLy*y)*sin(c_wxy + b_wxy*twopi_invLx*x) - a_wxz*b_wxz*twopi_invLx*cos(g_wxz + f_wxz*t)*cos(e_wxz + d_wxz*twopi_invLz*z)*sin(c_wxz + b_wxz*twopi_invLx*x))*cos(g_uyz + f_uyz*t)*sin(c_uyz + b_uyz*twopi_invLy*y)*sin(e_uyz + d_uyz*twopi_invLz*z) + 2*a_wxy*b_wxy*d_wxy*twopi_invLx*twopi_invLy*(-a_uz*b_uz*twopi_invLz*cos(g_uz + f_uz*t)*sin(c_uz + b_uz*twopi_invLz*z) - a_uxz*d_uxz*twopi_invLz*cos(g_uxz + f_uxz*t)*cos(c_uxz + b_uxz*twopi_invLx*x)*sin(e_uxz + d_uxz*twopi_invLz*z) - a_uyz*d_uyz*twopi_invLz*cos(c_uyz + b_uyz*twopi_invLy*y)*cos(g_uyz + f_uyz*t)*sin(e_uyz + d_uyz*twopi_invLz*z))*cos(g_wxy + f_wxy*t)*sin(c_wxy + b_wxy*twopi_invLx*x)*sin(e_wxy + d_wxy*twopi_invLy*y) + 2*a_wyz*b_wyz*d_wyz*twopi_invLy*twopi_invLz*(-a_wz*b_wz*twopi_invLz*cos(g_wz + f_wz*t)*sin(c_wz + b_wz*twopi_invLz*z) - a_wxz*d_wxz*twopi_invLz*cos(g_wxz + f_wxz*t)*cos(c_wxz + b_wxz*twopi_invLx*x)*sin(e_wxz + d_wxz*twopi_invLz*z) - a_wyz*d_wyz*twopi_invLz*cos(g_wyz + f_wyz*t)*cos(c_wyz + b_wyz*twopi_invLy*y)*sin(e_wyz + d_wyz*twopi_invLz*z))*cos(g_wyz + f_wyz*t)*sin(c_wyz + b_wyz*twopi_invLy*y)*sin(e_wyz + d_wyz*twopi_invLz*z)


A2=A2.subs(-Ux*b_ux**2*twopi_invLx**2 - Uxy*b_uxy**2*twopi_invLx**2 - Uxz*b_uxz**2*twopi_invLx**2,collect((-Ux*b_ux**2*twopi_invLx**2 - Uxy*b_uxy**2*twopi_invLx**2 - Uxz*b_uxz**2*twopi_invLx**2),twopi_invLx**2))
A2=A2.subs(-Uxy*d_uxy**2*twopi_invLy**2 - Uy*b_uy**2*twopi_invLy**2 - Uyz*b_uyz**2*twopi_invLy**2,collect((-Uxy*d_uxy**2*twopi_invLy**2 - Uy*b_uy**2*twopi_invLy**2 - Uyz*b_uyz**2*twopi_invLy**2),twopi_invLy**2))
A2=A2.subs(-Wxy*d_wxy**2*twopi_invLy**2 - Wy*b_wy**2*twopi_invLy**2 - Wyz*b_wyz**2*twopi_invLy**2,collect((-Wxy*d_wxy**2*twopi_invLy**2 - Wy*b_wy**2*twopi_invLy**2 - Wyz*b_wyz**2*twopi_invLy**2),twopi_invLy**2))
A2=A2.subs(-Vx*b_vx**2*twopi_invLx**2 - Vxy*b_vxy**2*twopi_invLx**2 - Vxz*b_vxz**2*twopi_invLx**2,collect(-Vx*b_vx**2*twopi_invLx**2 - Vxy*b_vxy**2*twopi_invLx**2 - Vxz*b_vxz**2*twopi_invLx**2,twopi_invLx**2))
A2=A2.subs(-Wxz*d_wxz**2*twopi_invLz**2 - Wyz*d_wyz**2*twopi_invLz**2 - Wz*b_wz**2*twopi_invLz**2,collect(-Wxz*d_wxz**2*twopi_invLz**2 - Wyz*d_wyz**2*twopi_invLz**2 - Wz*b_wz**2*twopi_invLz**2,twopi_invLz**2))
A2=A2.subs(-Uxz*d_uxz**2*twopi_invLz**2 - Uyz*d_uyz**2*twopi_invLz**2 - Uz*b_uz**2*twopi_invLz**2,collect(-Uxz*d_uxz**2*twopi_invLz**2 - Uyz*d_uyz**2*twopi_invLz**2 - Uz*b_uz**2*twopi_invLz**2,twopi_invLz**2))
A2=A2.subs(-Vxz*d_vxz**2*twopi_invLz**2 - Vyz*d_vyz**2*twopi_invLz**2 - Vz*b_vz**2*twopi_invLz**2,collect(-Vxz*d_vxz**2*twopi_invLz**2 - Vyz*d_vyz**2*twopi_invLz**2 - Vz*b_vz**2*twopi_invLz**2,twopi_invLz**2))
A2=A2.subs(-Wx*b_wx**2*twopi_invLx**2 - Wxy*b_wxy**2*twopi_invLx**2 - Wxz*b_wxz**2*twopi_invLx**2,collect(-Wx*b_wx**2*twopi_invLx**2 - Wxy*b_wxy**2*twopi_invLx**2 - Wxz*b_wxz**2*twopi_invLx**2,twopi_invLx**2))


#collect(A2,[twopi_invLx**2*(-Ux*b_ux**2 - Uxy*b_uxy**2 - Uxz*b_uxz**2), twopi_invLy**2*(-Uxy*d_uxy**2 - Uy*b_uy**2 - Uyz*b_uyz**2),twopi_invLy**2*(-Wxy*d_wxy**2 - Wy*b_wy**2 - Wyz*b_wyz**2),twopi_invLz**2*(-Wxz*d_wxz**2 - Wyz*d_wyz**2 - Wz*b_wz**2),twopi_invLx**2*(-Wx*b_wx**2 - Wxy*b_wxy**2 - Wxz*b_wxz**2), twopi_invLz**2*(-Uxz*d_uxz**2 - Uyz*d_uyz**2 - Uz*b_uz**2),twopi_invLz**2*(-Vxz*d_vxz**2 - Vyz*d_vyz**2 - Vz*b_vz**2),twopi_invLx**2*(-Vx*b_vx**2 - Vxy*b_vxy**2 - Vxz*b_vxz**2)])
A2=collect(A2,[(-Ux*b_ux**2 - Uxy*b_uxy**2 - Uxz*b_uxz**2),(-Uxy*d_uxy**2 - Uy*b_uy**2 - Uyz*b_uyz**2),(-Wxy*d_wxy**2 - Wy*b_wy**2 - Wyz*b_wyz**2),(-Wxz*d_wxz**2 - Wyz*d_wyz**2 - Wz*b_wz**2),(-Wx*b_wx**2 - Wxy*b_wxy**2 - Wxz*b_wxz**2), (-Uxz*d_uxz**2 - Uyz*d_uyz**2 - Uz*b_uz**2),(-Vxz*d_vxz**2 - Vyz*d_vyz**2 - Vz*b_vz**2),(-Vx*b_vx**2 - Vxy*b_vxy**2 - Vxz*b_vxz**2)])





Q2n=A1+A2



Q2n=Q2n.subs(U0  , a_u0*cos(g_u0 + f_u0*t));
Q2n=Q2n.subs(Ux , a_ux*cos(c_ux + b_ux*twopi_invLx*x)*cos(g_ux + f_ux*t));
Q2n=Q2n.subs(Uy , a_uy*cos(g_uy + f_uy*t)*cos(c_uy + b_uy*twopi_invLy*y)); 
Q2n=Q2n.subs(Uz , a_uz*cos(g_uz + f_uz*t)*cos(c_uz + b_uz*twopi_invLz*z));
Q2n=Q2n.subs(Uxy , a_uxy*cos(c_uxy + b_uxy*twopi_invLx*x)*cos(e_uxy + d_uxy*twopi_invLy*y)*cos(g_uxy + f_uxy*t));
Q2n=Q2n.subs(Uxz , a_uxz*cos(g_uxz + f_uxz*t)*cos(c_uxz + b_uxz*twopi_invLx*x)*cos(e_uxz + d_uxz*twopi_invLz*z));
Q2n=Q2n.subs(Uyz , a_uyz*cos(c_uyz + b_uyz*twopi_invLy*y)*cos(g_uyz + f_uyz*t)*cos(e_uyz + d_uyz*twopi_invLz*z));

Q2n=Q2n.subs(V0 , a_v0*cos(g_v0 + f_v0*t));
Q2n=Q2n.subs(Vx , a_vx*cos(c_vx + b_vx*twopi_invLx*x)*cos(g_vx + f_vx*t));
Q2n=Q2n.subs(Vy , a_vy*cos(c_vy + b_vy*twopi_invLy*y)*cos(g_vy + f_vy*t));
Q2n=Q2n.subs(Vz , a_vz*cos(g_vz + f_vz*t)*cos(c_vz + b_vz*twopi_invLz*z));
Q2n=Q2n.subs(Vxy , a_vxy*cos(e_vxy + d_vxy*twopi_invLy*y)*cos(g_vxy + f_vxy*t)*cos(c_vxy + b_vxy*twopi_invLx*x));
Q2n=Q2n.subs(Vxz , a_vxz*cos(e_vxz + d_vxz*twopi_invLz*z)*cos(g_vxz + f_vxz*t)*cos(c_vxz + b_vxz*twopi_invLx*x));
Q2n=Q2n.subs(Vyz , a_vyz*cos(g_vyz + f_vyz*t)*cos(c_vyz + b_vyz*twopi_invLy*y)*cos(e_vyz + d_vyz*twopi_invLz*z));

Q2n=Q2n.subs(W0, a_w0*cos(g_w0 + f_w0*t) );
Q2n=Q2n.subs(Wx, a_wx*cos(g_wx + f_wx*t)*cos(c_wx + b_wx*twopi_invLx*x));
Q2n=Q2n.subs(Wy, a_wy*cos(c_wy + b_wy*twopi_invLy*y)*cos(g_wy + f_wy*t));
Q2n=Q2n.subs(Wz, a_wz*cos(g_wz + f_wz*t)*cos(c_wz + b_wz*twopi_invLz*z));
Q2n=Q2n.subs(Wxy, a_wxy*cos(g_wxy + f_wxy*t)*cos(c_wxy + b_wxy*twopi_invLx*x)*cos(e_wxy + d_wxy*twopi_invLy*y) );
Q2n=Q2n.subs(Wxz, a_wxz*cos(g_wxz + f_wxz*t)*cos(c_wxz + b_wxz*twopi_invLx*x)*cos(e_wxz + d_wxz*twopi_invLz*z) );
Q2n=Q2n.subs(Wyz, a_wyz*cos(g_wyz + f_wyz*t)*cos(c_wyz + b_wyz*twopi_invLy*y)*cos(e_wyz + d_wyz*twopi_invLz*z) );




#Checking Q2 --------------------------------------------------------
Q2new=Q2n.subs(U,a_u0*cos(g_u0 + f_u0*t) + a_ux*cos(c_ux + b_ux*twopi_invLx*x)*cos(g_ux + f_ux*t) + a_uy*cos(g_uy + f_uy*t)*cos(c_uy + b_uy*twopi_invLy*y) + a_uz*cos(g_uz + f_uz*t)*cos(c_uz + b_uz*twopi_invLz*z) + a_uxy*cos(c_uxy + b_uxy*twopi_invLx*x)*cos(e_uxy + d_uxy*twopi_invLy*y)*cos(g_uxy + f_uxy*t) + a_uxz*cos(g_uxz + f_uxz*t)*cos(c_uxz + b_uxz*twopi_invLx*x)*cos(e_uxz + d_uxz*twopi_invLz*z) + a_uyz*cos(c_uyz + b_uyz*twopi_invLy*y)*cos(g_uyz + f_uyz*t)*cos(e_uyz + d_uyz*twopi_invLz*z)
);
Q2new=Q2new.subs(V,a_v0*cos(g_v0 + f_v0*t) + a_vx*cos(c_vx + b_vx*twopi_invLx*x)*cos(g_vx + f_vx*t) + a_vy*cos(c_vy + b_vy*twopi_invLy*y)*cos(g_vy + f_vy*t) + a_vz*cos(g_vz + f_vz*t)*cos(c_vz + b_vz*twopi_invLz*z) + a_vxy*cos(e_vxy + d_vxy*twopi_invLy*y)*cos(g_vxy + f_vxy*t)*cos(c_vxy + b_vxy*twopi_invLx*x) + a_vxz*cos(e_vxz + d_vxz*twopi_invLz*z)*cos(g_vxz + f_vxz*t)*cos(c_vxz + b_vxz*twopi_invLx*x) + a_vyz*cos(g_vyz + f_vyz*t)*cos(c_vyz + b_vyz*twopi_invLy*y)*cos(e_vyz + d_vyz*twopi_invLz*z)
);
Q2new=Q2new.subs(W,a_w0*cos(g_w0 + f_w0*t) + a_wx*cos(g_wx + f_wx*t)*cos(c_wx + b_wx*twopi_invLx*x) + a_wy*cos(c_wy + b_wy*twopi_invLy*y)*cos(g_wy + f_wy*t) + a_wz*cos(g_wz + f_wz*t)*cos(c_wz + b_wz*twopi_invLz*z) + a_wxy*cos(g_wxy + f_wxy*t)*cos(c_wxy + b_wxy*twopi_invLx*x)*cos(e_wxy + d_wxy*twopi_invLy*y) + a_wxz*cos(g_wxz + f_wxz*t)*cos(c_wxz + b_wxz*twopi_invLx*x)*cos(e_wxz + d_wxz*twopi_invLz*z) + a_wyz*cos(g_wyz + f_wyz*t)*cos(c_wyz + b_wyz*twopi_invLy*y)*cos(e_wyz + d_wyz*twopi_invLz*z)
);
Q2new=expand(Q2new)
Q2=expand(Q2)

Q2new==Q2 # true

#=====================================================================


Q_v=Q1n+Q2n+Q3n

U=(a_u0*cos(g_u0 + f_u0*t) + a_ux*cos(c_ux + b_ux*twopi_invLx*x)*cos(g_ux + f_ux*t) + a_uy*cos(g_uy + f_uy*t)*cos(c_uy + b_uy*twopi_invLy*y) + a_uz*cos(g_uz + f_uz*t)*cos(c_uz + b_uz*twopi_invLz*z) + a_uxy*cos(c_uxy + b_uxy*twopi_invLx*x)*cos(e_uxy + d_uxy*twopi_invLy*y)*cos(g_uxy + f_uxy*t) + a_uxz*cos(g_uxz + f_uxz*t)*cos(c_uxz + b_uxz*twopi_invLx*x)*cos(e_uxz + d_uxz*twopi_invLz*z) + a_uyz*cos(c_uyz + b_uyz*twopi_invLy*y)*cos(g_uyz + f_uyz*t)*cos(e_uyz + d_uyz*twopi_invLz*z)
);
V=(a_v0*cos(g_v0 + f_v0*t) + a_vx*cos(c_vx + b_vx*twopi_invLx*x)*cos(g_vx + f_vx*t) + a_vy*cos(c_vy + b_vy*twopi_invLy*y)*cos(g_vy + f_vy*t) + a_vz*cos(g_vz + f_vz*t)*cos(c_vz + b_vz*twopi_invLz*z) + a_vxy*cos(e_vxy + d_vxy*twopi_invLy*y)*cos(g_vxy + f_vxy*t)*cos(c_vxy + b_vxy*twopi_invLx*x) + a_vxz*cos(e_vxz + d_vxz*twopi_invLz*z)*cos(g_vxz + f_vxz*t)*cos(c_vxz + b_vxz*twopi_invLx*x) + a_vyz*cos(g_vyz + f_vyz*t)*cos(c_vyz + b_vyz*twopi_invLy*y)*cos(e_vyz + d_vyz*twopi_invLz*z)
);
W=(a_w0*cos(g_w0 + f_w0*t) + a_wx*cos(g_wx + f_wx*t)*cos(c_wx + b_wx*twopi_invLx*x) + a_wy*cos(c_wy + b_wy*twopi_invLy*y)*cos(g_wy + f_wy*t) + a_wz*cos(g_wz + f_wz*t)*cos(c_wz + b_wz*twopi_invLz*z) + a_wxy*cos(g_wxy + f_wxy*t)*cos(c_wxy + b_wxy*twopi_invLx*x)*cos(e_wxy + d_wxy*twopi_invLy*y) + a_wxz*cos(g_wxz + f_wxz*t)*cos(c_wxz + b_wxz*twopi_invLx*x)*cos(e_wxz + d_wxz*twopi_invLz*z) + a_wyz*cos(g_wyz + f_wyz*t)*cos(c_wyz + b_wyz*twopi_invLy*y)*cos(e_wyz + d_wyz*twopi_invLz*z)
);
from sympy.utilities.codegen import codegen,Routine
codegen((
	
	("U", U),
	("V", V),
	("W", W),
        ("Q_v",    Q_v   ),

        ), "C", "../C_code_sympy/incompressible_flow_source_Qv", header=True, to_files=True)




# Writing Q_g into a latex file ------------------------------------------
latexQ_v=latex(Q_v)
s=str(latexQ_v)
f=open('../latex/Qv.tex','w')
f.write(s)
f.close()



# Writing Q_v into a file in order to count characters ------------------------------------------
s=str(Q_v)
f=open('../C_code_sympy/SourceQv_after_factorization.dat','w')
f.write(s)
f.close()
