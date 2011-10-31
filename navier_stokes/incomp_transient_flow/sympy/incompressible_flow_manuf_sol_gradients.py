
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

hg = diff(H1,z) - diff(H3,x)


# Eq_g = diff(g,t) == hg + (1/Re) * Laplacian g ----------------------------------------------
diff(g, t)+ hg == (1/Re)* (diff(g,x,2)+ diff(g, y,2)+  diff(g,z,2))

Lo= diff(g,t)+ hg -(1/Re)* (diff(g,x,2)+ diff(g, y,2)+ diff(g,z,2))

L1=diff(g,t);
L2=+hg
L3=-(1/Re)* (diff(g,x,2)+ diff(g, y,2)+ diff(g,z,2))

R=Lo-(L1+L2+L3)




# Manufactured solutions ------------------------------------------------
# Solution parameters used in the form of the analytical solution
var(""" a_0  a_x  a_xy  a_xz  a_y  a_yz  a_z
             b_x  b_xy  b_xz  b_y  b_yz  b_z
             c_x  c_xy  c_xz  c_y  c_yz  c_z
                  d_xy  d_xz       d_yz
                  e_xy  e_xz       e_yz
        f_0  f_x  f_xy  f_xz  f_y  f_yz  f_z
        g_0  g_x  g_xy  g_xz  g_y  g_yz  g_z """, real=True)

# Explicitly keep (2 * pi / L) terms together as indivisible tokens
var('twopi_invLx  twopi_invLy  twopi_invLz', real=True)
var(""" a_u0,a_ux,a_uxy,a_uxz,a_uy,a_uyz,a_uz,b_ux,b_uxy,b_uxz,b_uy,b_uyz,b_uz, c_ux,c_uxy,c_uxz,c_uy,c_uyz,c_uz, d_uxy,d_uxz,d_uyz, e_uxy,e_uxz,e_uyz, f_u0,f_ux,f_uxy,f_uxz,f_uy,f_uyz,f_uz, g_u0,g_ux,g_uxy,g_uxz,g_uy,g_uyz,g_uz""",real=True)
var(""" a_v0,a_vx,a_vxy,a_vxz,a_vy,a_vyz,a_vz,b_vx,b_vxy,b_vxz,b_vy,b_vyz,b_vz, c_vx,c_vxy,c_vxz,c_vy,c_vyz,c_vz, d_vxy,d_vxz,d_vyz, e_vxy,e_vxz,e_vyz, f_v0,f_vx,f_vxy,f_vxz,f_vy,f_vyz,f_vz, g_v0,g_vx,g_vxy,g_vxz,g_vy,g_vyz,g_vz""",real=True);
var(""" a_w0,a_wx,a_wxy,a_wxz,a_wy,a_wyz,a_wz,b_wx,b_wxy,b_wxz,b_wy,b_wyz,b_wz, c_wx,c_wxy,c_wxz,c_wy,c_wyz,c_wz, d_wxy,d_wxz,d_wyz, e_wxy,e_wxz,e_wyz, f_w0,f_wx,f_wxy,f_wxz,f_wy,f_wyz,f_wz, g_w0,g_wx,g_wxy,g_wxz,g_wy,g_wyz,g_wz""",real=True)

# Form the analytical solution and its derivatives
phi  = (
        a_0                                                                 *cos(f_0 *t + g_0 )
      + a_x  * cos(b_x *twopi_invLx*x + c_x )                               *cos(f_x *t + g_x )
      + a_xy * cos(b_xy*twopi_invLx*x + c_xy)*cos(d_xy*twopi_invLy*y + e_xy)*cos(f_xy*t + g_xy)
      + a_xz * cos(b_xz*twopi_invLx*x + c_xz)*cos(d_xz*twopi_invLz*z + e_xz)*cos(f_xz*t + g_xz)
      + a_y  * cos(b_y *twopi_invLy*y + c_y )                               *cos(f_y *t + g_y )
      + a_yz * cos(b_yz*twopi_invLy*y + c_yz)*cos(d_yz*twopi_invLz*z + e_yz)*cos(f_yz*t + g_yz)
      + a_z  * cos(b_z *twopi_invLz*z + c_z )                               *cos(f_z *t + g_z )
)

params_u = {   
    'a_0'  : a_u0,    'a_x'  : a_ux,     'a_xy' : a_uxy,    'a_xz' : a_uxz,
    'a_y'  : a_uy,    'a_yz' : a_uyz,    'a_z'  : a_uz,
    'b_x'  : b_ux,    'b_xy' : b_uxy,    'b_xz' : b_uxz,
    'b_y'  : b_uy,    'b_yz' : b_uyz,    'b_z'  : b_uz,
    'c_x'  : c_ux,    'c_xy' : c_uxy,    'c_xz' : c_uxz,
    'c_y'  : c_uy,    'c_yz' : c_uyz,    'c_z'  : c_uz,
    'd_xy' : d_uxy,   'd_xz' : d_uxz,    'd_yz' : d_uyz,
    'e_xy' : e_uxy,   'e_xz' : e_uxz,    'e_yz' : e_uyz,
    'f_0'  : f_u0,    'f_x'  : f_ux,     'f_xy' : f_uxy,
    'f_xz' : f_uxz,   'f_y'  : f_uy,     'f_yz' : f_uyz,
    'f_z'  : f_uz,    'g_0'  : g_u0,     'g_x'  : g_ux,
    'g_xy' : g_uxy,   'g_xz' : g_uxz,    'g_y'  : g_uy,
    'g_yz' : g_uyz,   'g_z'  : g_uz,
}
params_v = {   
    'a_0'  : a_v0,    'a_x'  : a_vx,     'a_xy' : a_vxy,    'a_xz' : a_vxz,
    'a_y'  : a_vy,    'a_yz' : a_vyz,    'a_z'  : a_vz,
    'b_x'  : b_vx,    'b_xy' : b_vxy,    'b_xz' : b_vxz,
    'b_y'  : b_vy,    'b_yz' : b_vyz,    'b_z'  : b_vz,
    'c_x'  : c_vx,    'c_xy' : c_vxy,    'c_xz' : c_vxz,
    'c_y'  : c_vy,    'c_yz' : c_vyz,    'c_z'  : c_vz,
    'd_xy' : d_vxy,   'd_xz' : d_vxz,    'd_yz' : d_vyz,
    'e_xy' : e_vxy,   'e_xz' : e_vxz,    'e_yz' : e_vyz,
    'f_0'  : f_v0,    'f_x'  : f_vx,     'f_xy' : f_vxy,
    'f_xz' : f_vxz,   'f_y'  : f_vy,     'f_yz' : f_vyz,
    'f_z'  : f_vz,    'g_0'  : g_v0,     'g_x'  : g_vx,
    'g_xy' : g_vxy,   'g_xz' : g_vxz,    'g_y'  : g_vy,
    'g_yz' : g_vyz,   'g_z'  : g_vz,
}
params_w = {
    'a_0'  : a_w0,    'a_x'  : a_wx,     'a_xy' : a_wxy,    'a_xz' : a_wxz,
    'a_y'  : a_wy,    'a_yz' : a_wyz,    'a_z'  : a_wz,
    'b_x'  : b_wx,    'b_xy' : b_wxy,    'b_xz' : b_wxz,
    'b_y'  : b_wy,    'b_yz' : b_wyz,    'b_z'  : b_wz,
    'c_x'  : c_wx,    'c_xy' : c_wxy,    'c_xz' : c_wxz,
    'c_y'  : c_wy,    'c_yz' : c_wyz,    'c_z'  : c_wz,
    'd_xy' : d_wxy,   'd_xz' : d_wxz,    'd_yz' : d_wyz,
    'e_xy' : e_wxy,   'e_xz' : e_wxz,    'e_yz' : e_wyz,
    'f_0'  : f_w0,    'f_x'  : f_wx,     'f_xy' : f_wxy,
    'f_xz' : f_wxz,   'f_y'  : f_wy,     'f_yz' : f_wyz,
    'f_z'  : f_wz,    'g_0'  : g_w0,     'g_x'  : g_wx,
    'g_xy' : g_wxy,   'g_xz' : g_wxz,    'g_y'  : g_wy,
    'g_yz' : g_wyz,   'g_z'  : g_wz,
}

u_an  = phi .subs(params_u)
v_an  = phi .subs(params_v)
w_an  = phi .subs(params_w)

g_an = diff(u_an,z)-diff(w_an,x)

grad_u_an=[diff(u_an,x),diff(u_an,y),diff(u_an,z)]

# Writing u_an,v_an, and w_an  and their gradients in a C file ------------------------------------------
dudx=diff(u_an,x);
dudy=diff(u_an,y);
dudz=diff(u_an,z);


dvdx=diff(v_an,x);
dvdy=diff(v_an,y);
dvdz=diff(v_an,z);


dwdx=diff(w_an,x);
dwdy=diff(w_an,y);
dwdz=diff(w_an,z);


from sympy.printing.ccode import ccode, CCodePrinter
ccode(u_an)



from sympy.utilities.codegen import codegen,Routine
codegen((
	
	("u_an", u_an),
	("v_an", v_an),
	("w_an", w_an),
	("du_dx",dudx),
	("du_dy",dudy),
	("du_dz",dudz),
	("dv_dx",dvdx),
	("dv_dy",dvdy),
	("dv_dz",dvdz),
        ("dw_dx",dwdx),
	("dw_dy",dwdy),
	("dw_dz",dwdz),
        
        ), "C", "../C_code_sympy/incompressible_flow_manuf_solutions_gradients", header=True, to_files=True)



# Writing Q_g into a latex file ------------------------------------------
f=open('../latex/Derivatives.tex','w')
f.write(latex(dudx));
f.write(latex(dudy));
f.write(latex(dudz));
f.write(latex(dvdx));
f.write(latex(dvdy));
f.write(latex(dvdz));
f.write(latex(dwdx));
f.write(latex(dwdy));
f.write(latex(dwdz));


f.close()



#import numpy as np 
#gradu = np.array([diff(u_an,x),diff(u_an,y),diff(u_an,z)])

#def Sourceteste():
  #U=a_u0*cos(g_u0 + f_u0*t)*cos(g_u0 + f_u0*t)
  #V=a_v0*cos(g_v0 + f_v0*t)
  #W=U+V 
  #return W

#w= Sourceteste()
##In [59]:  np.array([1, 2, 3])
