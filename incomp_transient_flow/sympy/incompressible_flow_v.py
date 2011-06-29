
from sympy import *

var('x y z t')

var('Re')

u = Function('u')(x,y,z,t)
v = Function('v')(x,y,z,t)
w = Function('w')(x,y,z,t)
g = Function('g')(x,y,z,t)

hv =Function('hv')(x,y,z,t)
hg =Function('hg')(x,y,z,t)
# =Function(' ')(x,y,z,t)

H1 = u*(diff(u, x))+v*(diff(u, y))+w*(diff(u, z))
H2 = u*(diff(v, x))+v*(diff(v, y))+w*(diff(v, z))
H3 = u*(diff(w, x))+v*(diff(w, y))+w*(diff(w, z))

aux1= diff(H1,x) + diff(H3,z)

hv=-diff(aux1,y) + diff(H2,x,x,)+ diff(H2,z,z)

# Eq_v = diff(Laplacian v,t) == hv+1/Re (Laplacian^2 v) -----------------------------------------------------
diff(diff(v, x, x)+ diff(v, y, y)+ diff(v, z, z), t)== hv+(1/Re)* (diff(v,x,4)+ 2*diff(v,y,y,x,x)+ 2*diff(v,z,z,x,x)+ diff(v, y,4)+ 2*diff(v, z, z, y, y)+ diff(v,z,4))

Lo=diff(diff(v, x, x)+ diff(v, y, y)+ diff(v, z, z), t)- hv-(1/Re)* (diff(v,x,4)+ 2*diff(v,y,y,x,x)+ 2*diff(v,z,z,x,x)+ diff(v, y,4)+ 2*diff(v, z, z, y, y)+ diff(v,z,4))

L1=diff(diff(v, x, x)+ diff(v, y, y)+ diff(v, z, z), t);
L2=-hv
L3=-(1/Re)* (diff(v,x,4)+ 2*diff(v,y,y,x,x)+ 2*diff(v,z,z,x,x)+ diff(v, y,4)+ 2*diff(v, z, z, y, y)+ diff(v,z,4))

R=Lo-(L1+L2+L3)


# Manufactured solutions ----------------------------------------------------------------------------------------

# Explicitly keep (2 * pi / L) terms together as indivisible tokens
var('twopi_invLx  twopi_invLy  twopi_invLz', real=True)
var(""" a_u0,a_ux,a_uxy,a_uxz,a_uy,a_uyz,a_uz,b_ux,b_uxy,b_uxz,b_uy,b_uyz,b_uz, c_ux,c_uxy,c_uxz,c_uy,c_uyz,c_uz, d_uxy,d_uxz,d_uyz, e_uxy,e_uxz,e_uyz, f_u0,f_ux,f_uxy,f_uxz,f_uy,f_uyz,f_uz, g_u0,g_ux,g_uxy,g_uxz,g_uy,g_uyz,g_uz""",real=True)
var(""" a_v0,a_vx,a_vxy,a_vxz,a_vy,a_vyz,a_vz,b_vx,b_vxy,b_vxz,b_vy,b_vyz,b_vz, c_vx,c_vxy,c_vxz,c_vy,c_vyz,c_vz, d_vxy,d_vxz,d_vyz, e_vxy,e_vxz,e_vyz, f_v0,f_vx,f_vxy,f_vxz,f_vy,f_vyz,f_vz, g_v0,g_vx,g_vxy,g_vxz,g_vy,g_vyz,g_vz""",real=True);
var(""" a_w0,a_wx,a_wxy,a_wxz,a_wy,a_wyz,a_wz,b_wx,b_wxy,b_wxz,b_wy,b_wyz,b_wz, c_wx,c_wxy,c_wxz,c_wy,c_wyz,c_wz, d_wxy,d_wxz,d_wyz, e_wxy,e_wxz,e_wyz, f_w0,f_wx,f_wxy,f_wxz,f_wy,f_wyz,f_wz, g_w0,g_wx,g_wxy,g_wxz,g_wy,g_wyz,g_wz""",real=True)

var("u_an, v_an, w_an")

u_an = a_u0*cos(g_u0 + f_u0*t) + a_ux*cos(c_ux + b_ux*twopi_invLx*x)*cos(g_ux + f_ux*t) + a_uy*cos(g_uy + f_uy*t)*cos(c_uy + b_uy*twopi_invLy*y) + a_uz*cos(g_uz + f_uz*t)*cos(c_uz + b_uz*twopi_invLz*z) + a_uxy*cos(c_uxy + b_uxy*twopi_invLx*x)*cos(e_uxy + d_uxy*twopi_invLy*y)*cos(g_uxy + f_uxy*t) + a_uxz*cos(g_uxz + f_uxz*t)*cos(c_uxz + b_uxz*twopi_invLx*x)*cos(e_uxz + d_uxz*twopi_invLz*z) + a_uyz*cos(c_uyz + b_uyz*twopi_invLy*y)*cos(g_uyz + f_uyz*t)*cos(e_uyz + d_uyz*twopi_invLz*z);
v_an = a_v0*cos(g_v0 + f_v0*t) + a_vx*cos(c_vx + b_vx*twopi_invLx*x)*cos(g_vx + f_vx*t) + a_vy*cos(c_vy + b_vy*twopi_invLy*y)*cos(g_vy + f_vy*t) + a_vz*cos(g_vz + f_vz*t)*cos(c_vz + b_vz*twopi_invLz*z) + a_vxy*cos(e_vxy + d_vxy*twopi_invLy*y)*cos(g_vxy + f_vxy*t)*cos(c_vxy + b_vxy*twopi_invLx*x) + a_vxz*cos(e_vxz + d_vxz*twopi_invLz*z)*cos(g_vxz + f_vxz*t)*cos(c_vxz + b_vxz*twopi_invLx*x) + a_vyz*cos(g_vyz + f_vyz*t)*cos(c_vyz + b_vyz*twopi_invLy*y)*cos(e_vyz + d_vyz*twopi_invLz*z);
w_an = a_w0*cos(g_w0 + f_w0*t) + a_wx*cos(g_wx + f_wx*t)*cos(c_wx + b_wx*twopi_invLx*x) + a_wy*cos(c_wy + b_wy*twopi_invLy*y)*cos(g_wy + f_wy*t) + a_wz*cos(g_wz
+ f_wz*t)*cos(c_wz + b_wz*twopi_invLz*z) + a_wxy*cos(g_wxy + f_wxy*t)*cos(c_wxy + b_wxy*twopi_invLx*x)*cos(e_wxy + d_wxy*twopi_invLy*y) + a_wxz*cos(g_wxz +
f_wxz*t)*cos(c_wxz + b_wxz*twopi_invLx*x)*cos(e_wxz + d_wxz*twopi_invLz*z) + a_wyz*cos(g_wyz + f_wyz*t)*cos(c_wyz + b_wyz*twopi_invLy*y)*cos(e_wyz +
d_wyz*twopi_invLz*z);


# Aplying L1 on v in order to obtain source term Q1 ----------

Q1=L1.subs(v,v_an)

# Aplying L2 on u, v and w in order to obtain source term Q2 ----------
Q2=L2.subs(u,u_an)
Q2=Q2.subs(v,v_an)
Q2=Q2.subs(w,w_an)

# Aplying L3 on u, v and w in order to obtain source term Q3 ----------
Q3=L3.subs(u,u_an)
Q3=Q3.subs(v,v_an)
Q3=Q3.subs(w,w_an)


Qv=Q1+Q2+Q3



# Writing Q_v into a file in order to count characters ------------------------------------------
s=str(Qv)
f=open('../C_code_sympy/SourceQv_before_factorization.dat','w')
f.write(s)
f.close()
