# This program numerically checks that the manipulations done over source term Q for the 3D transient Navier-Stokes equations with Power Law
# viscosity model are correct. 
# mu = A_mu*T^(3/2)/(T+B_mu).

# Velocity u (momentum)

from sympy import *

var('x y z t')

u = Function('u')(x,y,z,t)
v = Function('v')(x,y,z,t)
w = Function('w')(x,y,z,t)
rho = Function('rho')(x,y,z,t)
p = Function('p')(x,y,z,t)
T = Function('T')(x,y,z,t)
e_t = Function('e_t')(x,y,z,t)


mu = Function('mu')(x,y,z,t)

tau_xx = Function('tau_xx')(x,y,z,t)
tau_xy = Function('tau_xy')(x,y,z,t)
tau_xz = Function('tau_xz')(x,y,z,t)
tau_yx = Function('tau_yx')(x,y,z,t)
tau_yy = Function('tau_yy')(x,y,z,t)
tau_yz = Function('tau_yz')(x,y,z,t)
tau_zx = Function('tau_zx')(x,y,z,t)
tau_zy = Function('tau_zy')(x,y,z,t)
tau_zz = Function('tau_zz')(x,y,z,t)


#3D transient Navier-Stokes equation - velocity u:
diff(rho*u,t) + diff(rho*u**2,x) + diff(rho*u*v,y) + diff(rho*u*w,z) + diff(p,x) ==diff(tau_xx, x)+ diff(tau_xy,y) + diff(tau_xz,z);
  


# Auxiliary relations -------------------------------------------------------------------------------------------
var("R,k,gamma,E_t,mu_ref,T_ref,beta,alpha")
T = p/(rho*R);
e = R*T/(gamma-1);
e_t = e+(u*u+v*v+w*w)/2;


mu = mu_ref * (T/T_ref)**beta; #it was commented in the other files
Lambda = (alpha- Rational(2,3) )*mu;

tau_xx = 2 * mu * diff(u, x) + Lambda * ( diff(u, x) + diff(v, y) + diff(w, z) );
tau_yy = 2 * mu * diff(v, y) + Lambda * ( diff(u, x) + diff(v, y) + diff(w, z) );
tau_zz = 2 * mu * diff(w, z) + Lambda * ( diff(u, x) + diff(v, y) + diff(w, z) );
tau_xy = mu*(diff(u, y)+diff(v, x));
tau_xz = mu*(diff(u, z)+diff(w, x));
tau_yz = mu*(diff(w, y)+diff(v, z));
tau_yx = tau_xy;
tau_zx = tau_xz;
tau_zy = tau_yz;

# Velocity-u equation, writen as differential operator:
Lo =diff(rho*u,t) + diff(rho*u**2,x) + diff(rho*u*v,y) + diff(rho*u*w,z) + diff(p,x) -(diff(tau_xx, x)+ diff(tau_xy,y) + diff(tau_xz,z));

L1 = +diff(rho*u**2, x)+diff(rho*u*v, y)+diff(rho*u*w, z);
L2 = +diff(p, x);
L3 = -(diff(tau_xx, x) + diff(tau_xy,y) + diff(tau_xz,z));
L4 = +diff(rho*u, t);

Res=Lo-L1-L2-L3-L4;
Res=expand(Res)
Res== 0 #True
#print(Res)

# Manufactured solutions ----------------------------------------------------------------------------------------
# Manufactured solutions ----------------------------------------------------------------------------------------
var("""L,Lt,rho_0,rho_x,a_rhox,rho_y,a_rhoy,rho_z,a_rhoz,rho_t,a_rhot, u_0,u_x,a_ux,u_y,a_uy,u_z,a_uz,u_t,a_ut,
v_0,v_x,a_vx,v_y,a_vy,v_z,a_vz,v_t,a_vt,w_0,w_x,a_wx,w_y,a_wy,w_z,a_wz,w_t,a_wt,p_0,p_x,a_px,p_y,a_py,p_z,a_pz,p_t,a_pt,
phi_0,phi_x,a_phix,phi_y,a_phiy,phi_z,a_phiz,phi_t,a_phit,""",real=True)



rho_an = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi * t
/ Lt);
u_an = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt);
v_an = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt);
w_an = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt);
p_an = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt);


# Convection: Aplying L1 on v in order to obtain source term Q1 ----------
Q=Lo.subs(u,u_an)
Q=Q.subs(v,v_an)
Q=Q.subs(w,w_an)
Q=Q.subs(rho,rho_an)
Qu=Q.subs(p,p_an)



# Work in (possibly slow) extended precision
PREC = 40;

# Use prime numbers for test data.  Slick generator for primes taken from
# http://stackoverflow.com/questions/567222/simple-prime-generator-in-python
def prime_generator():
    D = {}
    q = 2
    while True:
        if q not in D:
            yield q
            D[q * q] = [q]
        else:
            for p in D[q]:
                D.setdefault(p + q, []).append(p)
            del D[q]
        q += 1

# Test solution parameters.  These are not physically
# realizable but can be used to check the implementation.
p = prime_generator()
parameters = {
     gamma     : Real(p.next(), PREC),
     R         : Real(p.next(), PREC),
     mu_ref    : Real(p.next(), PREC),
     T_ref     : Real(p.next(), PREC),
     alpha     : Real(p.next(), PREC),
     beta     : Real(p.next(), PREC),
     k 	       : Real(p.next(), PREC),
     a_rhox    : Real(p.next(), PREC),
     a_rhoy    : Real(p.next(), PREC),
     a_rhoz    : Real(p.next(), PREC),
     a_rhot    : Real(p.next(), PREC),
     a_phix    : Real(p.next(), PREC),
     a_phiy    : Real(p.next(), PREC),
     a_phiz    : Real(p.next(), PREC),
     a_phit    : Real(p.next(), PREC),
    
     a_ux    : Real(p.next(), PREC),
     a_uy    : Real(p.next(), PREC),
     a_uz    : Real(p.next(), PREC),
     a_ut    : Real(p.next(), PREC),
    
     a_vx    : Real(p.next(), PREC),
     a_vy    : Real(p.next(), PREC),
     a_vz    : Real(p.next(), PREC),
     a_vt    : Real(p.next(), PREC),
    
     a_wx    : Real(p.next(), PREC),
     a_wy    : Real(p.next(), PREC),
     a_wz    : Real(p.next(), PREC),
     a_wt    : Real(p.next(), PREC),
    
     a_px    : Real(p.next(), PREC),
     a_py    : Real(p.next(), PREC),
     a_pz    : Real(p.next(), PREC),
     a_pt    : Real(p.next(), PREC),
     rho_0    : Real(p.next(), PREC),
     rho_x    : Real(p.next(), PREC),
     rho_y    : Real(p.next(), PREC),
     rho_z    : Real(p.next(), PREC),
     rho_t    : Real(p.next(), PREC),
     phi_0    : Real(p.next(), PREC),
     phi_x    : Real(p.next(), PREC),
     phi_y    : Real(p.next(), PREC),
     phi_z    : Real(p.next(), PREC),
     phi_t    : Real(p.next(), PREC),
     u_0    : Real(p.next(), PREC),
     u_x    : Real(p.next(), PREC),
     u_y    : Real(p.next(), PREC),
     u_z    : Real(p.next(), PREC),
     u_t    : Real(p.next(), PREC),
     v_0    : Real(p.next(), PREC),
     v_x    : Real(p.next(), PREC),
     v_y    : Real(p.next(), PREC),
     v_z    : Real(p.next(), PREC),
     v_t    : Real(p.next(), PREC),
     w_0    : Real(p.next(), PREC),
     w_x    : Real(p.next(), PREC),
     w_y    : Real(p.next(), PREC),
     w_z    : Real(p.next(), PREC),
     w_t    : Real(p.next(), PREC),
     p_0    : Real(p.next(), PREC),
     p_x    : Real(p.next(), PREC),
     p_y    : Real(p.next(), PREC),
     p_z    : Real(p.next(), PREC),
     p_t    : Real(p.next(), PREC),
     L        : Real(p.next(), PREC),
     Lt        : Real(p.next(), PREC),
     x : Rational(5, 10),
     y : Rational(6, 10),
     z : Rational(8, 10),
     t : Rational(9, 10),
    }

Qu_num=Qu.subs(parameters).evalf()  

Rho=rho_an;
U=u_an;
V=v_an;
W=w_an;
P=p_an;
T= P/(Rho*R);
Mu = mu_ref * (T/T_ref)**beta;

DMu_Dx=diff(Mu,x)
DMu_Dy=diff(Mu,y)
DMu_Dz=diff(Mu,z)

# it works either way. Note: needs to  either add a \ at the end, or put the command between parenthesis  
#DMu_Dx= (A_mu*T**0.5*(-1.5*pi*T*a_rhox*rho_x*cos(pi*a_rhox*x/L)/(L*Rho) - 1.5*pi*a_px*p_x*sin(pi*a_px*x/L)/(L*R*Rho))/(B_mu + T)
#+A_mu*T**1.5*(pi*T*a_rhox*rho_x*cos(pi*a_rhox*x/L)/(L*Rho) + pi*a_px*p_x*sin(pi*a_px*x/L)/(L*R*Rho))/(B_mu + T)**2 );

#DMu_Dy=( A_mu*T**0.5*(1.5*pi*T*a_rhoy*rho_y*sin(pi*a_rhoy*y/L)/(L*Rho) + 1.5*pi*a_py*p_y*cos(pi*a_py*y/L)/(L*R*Rho))/(B_mu + T)
#+A_mu*T**1.5*(-pi*T*a_rhoy*rho_y*sin(pi*a_rhoy*y/L)/(L*Rho) - pi*a_py*p_y*cos(pi*a_py*y/L)/(L*R*Rho))/(B_mu + T)**2 );

#DMu_Dz=( A_mu*T**0.5*(-1.5*pi*T*a_rhoz*rho_z*cos(pi*a_rhoz*z/L)/(L*Rho) - 1.5*pi*a_pz*p_z*sin(pi*a_pz*z/L)/(L*R*Rho))/(B_mu + T)
#+A_mu*T**1.5*(pi*T*a_rhoz*rho_z*cos(pi*a_rhoz*z/L)/(L*Rho) + pi*a_pz*p_z*sin(pi*a_pz*z/L)/(L*R*Rho))/(B_mu + T)**2 );

DMu_Dx=DMu_Dx.subs(parameters).evalf()  
DMu_Dy=DMu_Dy.subs(parameters).evalf()  
DMu_Dz=DMu_Dz.subs(parameters).evalf()

# continuation in line: either add a \ at the end, or put the command between parenthesis    
Q_u= ( pi*DMu_Dx*(-4*a_ux*u_x*cos(pi*a_ux*x/L)/3 - 2*a_wz*w_z*sin(pi*a_wz*z/L)/3 + 2*a_vy*v_y*cos(pi*a_vy*y/L)/3)/L +
pi*DMu_Dy*(a_uy*u_y*sin(pi*a_uy*y/L) + a_vx*v_x*sin(pi*a_vx*x/L))/L + pi*DMu_Dz*(a_uz*u_z*sin(pi*a_uz*z/L) - a_wx*w_x*cos(pi*a_wx*x/L))/L +
Mu*pi**2*(u_y*a_uy**2*cos(pi*a_uy*y/L) + u_z*a_uz**2*cos(pi*a_uz*z/L) + 4*u_x*a_ux**2*sin(pi*a_ux*x/L)/3)/L**2 + pi*Rho*U*(a_vy*v_y*cos(pi*a_vy*y/L)
-a_wz*w_z*sin(pi*a_wz*z/L) + 2*a_ux*u_x*cos(pi*a_ux*x/L))/L - pi*a_px*p_x*sin(pi*a_px*x/L)/L + pi*U*a_rhot*rho_t*cos(pi*a_rhot*t/Lt)/Lt +
pi*a_rhox*rho_x*U**2*cos(pi*a_rhox*x/L)/L - pi*Rho*a_ut*u_t*sin(pi*a_ut*t/Lt)/Lt + pi*U*W*a_rhoz*rho_z*cos(pi*a_rhoz*z/L)/L -
pi*Rho*V*a_uy*u_y*sin(pi*a_uy*y/L)/L - pi*Rho*W*a_uz*u_z*sin(pi*a_uz*z/L)/L - pi*U*V*a_rhoy*rho_y*sin(pi*a_rhoy*y/L)/L );


Q_u_num=Q_u.subs(parameters).evalf()  
  
  
print 'Q_u_num-Qu_num=',Q_u_num-Qu_num 
