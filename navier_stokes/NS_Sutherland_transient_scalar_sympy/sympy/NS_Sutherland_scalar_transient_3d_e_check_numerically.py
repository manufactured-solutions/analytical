# This program numerically checks that the manipulations done over source term Q for the 3D transient Navier-Stokes equations with Sutherland
# viscosity model are correct. 
# mu = A_mu*T^(3/2)/(T+B_mu).


# Total energy et

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

q_x = Function('q_x')(x,y,z,t)#(T)
q_y = Function('q_y')(x,y,z,t)#(T)
q_z = Function('q_z')(x,y,z,t)#(T)


#3D transient total energy equation
diff(rho*e_t, t)+diff(rho*u*e_t+p*u-u*tau_xx-v*tau_xy-w*tau_xz+q_x, x)+diff(rho*v*e_t+p*v-u*tau_yx-v*tau_yy-w*tau_yz+q_y,
y)+diff(rho*w*e_t+p*w-u*tau_zx-v*tau_zy-w*tau_zz+q_z, z)==0;



# Auxiliary relations -------------------------------------------------------------------------------------------
var("R,A_mu,B_mu,k,gamma,T_ref,mu_ref,Pr")
T = p/(rho*R);
e = R*T/(gamma-1);
e_t = e+(u*u+v*v+w*w)/2;
mu = A_mu*T**(3/2)/(T+B_mu);
k = gamma*R*mu/(gamma-1)/Pr;
tau_xx = mu*( Rational(4,3)*diff(u, x) -Rational(2,3)*diff(v, y) -Rational(2,3)*diff(w, z));
tau_yy = mu*( Rational(4,3)*diff(v, y) -Rational(2,3)*diff(u, x) -Rational(2,3)*diff(w, z));
tau_zz = mu*( Rational(4,3)*diff(w, z) -Rational(2,3)*diff(u, x) -Rational(2,3)*diff(v, y));
tau_xy = mu*(diff(u, y)+diff(v, x));
tau_xz = mu*(diff(u, z)+diff(w, x));
tau_yz = mu*(diff(w, y)+diff(v, z));
tau_yx = tau_xy;
tau_zx = tau_xz;
tau_zy = tau_yz;
q_x = -k*(diff(T, x));
q_y = -k*(diff(T, y));
q_z = -k*(diff(T, z));

# Energy equation, writen as differential operator:
Lo=diff(rho*e_t, t)+diff(rho*u*e_t+p*u-u*tau_xx-v*tau_xy-w*tau_xz+q_x, x)+diff(rho*v*e_t+p*v-u*tau_yx-v*tau_yy-w*tau_yz+q_y,
y)+diff(rho*w*e_t+p*w-u*tau_zx-v*tau_zy-w*tau_zz+q_z, z);  

L1=diff(rho*u*e_t,x) + diff(rho*v*e_t,y)+ diff(rho*w*e_t,z);
L2=+diff(p*u,x)+diff(p*v,y)+diff(p*w,z);
L3=+diff(-u*tau_xx-v*tau_xy-w*tau_xz,x)+diff(-u*tau_yx-v*tau_yy-w*tau_yz,y)+diff(-u*tau_zx-v*tau_zy-w*tau_zz,z);
L4=+diff(q_x, x)+diff(q_y, y)+diff(q_z, z);
L5=diff(rho*e_t, t);
 
 
Res=Lo-(L1+L2+L3+L4+L5)
Res=expand(Res)
Res== 0 #True

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


# Original source term: Aplying Lo on the manufactured solutions in order to obtain source term Qet ----------
Q=Lo.subs(u,u_an)
Q=Q.subs(v,v_an)
Q=Q.subs(w,w_an)
Q=Q.subs(rho,rho_an)
Qet=Q.subs(p,p_an)



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
     A_mu      : Real(p.next(), PREC),
     B_mu      : Real(p.next(), PREC),
     Pr        : Real(p.next(), PREC),
     
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

Qet_num=Qet.subs(parameters).evalf()  



Rho=rho_an;
U=u_an;
V=v_an;
W=w_an;
P=p_an;
E_t=-P/(Rho*(1 - gamma)) + U**2/2 + V**2/2 + W**2/2;
T= P/(Rho*R);

Mu = A_mu*T**(3/2)/(T+B_mu);

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

kappa = gamma*R*Mu/(gamma-1)/Pr;
kappa=kappa.subs(parameters).evalf()
# continuation in line: either add a \ at the end, or put the command between parenthesis    

Q_et=( pi*P*(a_ux*u_x*cos(pi*a_ux*x/L) + a_vy*v_y*cos(pi*a_vy*y/L) - a_wz*w_z*sin(pi*a_wz*z/L))/L + pi*Rho*(a_ux*u_x*U**2*cos(pi*a_ux*x/L) +
a_vy*v_y*V**2*cos(pi*a_vy*y/L) - a_wz*w_z*W**2*sin(pi*a_wz*z/L) + U*W*a_wx*w_x*cos(pi*a_wx*x/L) + V*W*a_vz*v_z*cos(pi*a_vz*z/L) +
V*W*a_wy*w_y*cos(pi*a_wy*y/L) - U*V*a_uy*u_y*sin(pi*a_uy*y/L) - U*V*a_vx*v_x*sin(pi*a_vx*x/L) - U*W*a_uz*u_z*sin(pi*a_uz*z/L))/L +
pi*Rho*(V*a_vt*v_t*cos(pi*a_vt*t/Lt) - U*a_ut*u_t*sin(pi*a_ut*t/Lt) - W*a_wt*w_t*sin(pi*a_wt*t/Lt))/Lt + pi*(U*a_px*p_x*sin(pi*a_px*x/L) +
W*a_pz*p_z*sin(pi*a_pz*z/L) - V*a_py*p_y*cos(pi*a_py*y/L))/(L*(1 - gamma)) + Mu*pi**2*(U*(u_y*a_uy**2*cos(pi*a_uy*y/L) + u_z*a_uz**2*cos(pi*a_uz*z/L)
+ 4*u_x*a_ux**2*sin(pi*a_ux*x/L)/3) + V*(v_x*a_vx**2*cos(pi*a_vx*x/L) + v_z*a_vz**2*sin(pi*a_vz*z/L) + 4*v_y*a_vy**2*sin(pi*a_vy*y/L)/3) +
W*(w_x*a_wx**2*sin(pi*a_wx*x/L) + w_y*a_wy**2*sin(pi*a_wy*y/L) + 4*w_z*a_wz**2*cos(pi*a_wz*z/L)/3) -
2*a_uy*a_vx*u_y*v_x*sin(pi*a_uy*y/L)*sin(pi*a_vx*x/L) - 2*a_vz*a_wy*v_z*w_y*cos(pi*a_vz*z/L)*cos(pi*a_wy*y/L) +
2*a_uz*a_wx*u_z*w_x*cos(pi*a_wx*x/L)*sin(pi*a_uz*z/L) - 4*a_ux*a_wz*u_x*w_z*cos(pi*a_ux*x/L)*sin(pi*a_wz*z/L)/3 -
4*a_vy*a_wz*v_y*w_z*cos(pi*a_vy*y/L)*sin(pi*a_wz*z/L)/3 + 4*a_ux*a_vy*u_x*v_y*cos(pi*a_ux*x/L)*cos(pi*a_vy*y/L)/3 - a_uy**2*u_y**2*sin(pi*a_uy*y/L)**2
- a_uz**2*u_z**2*sin(pi*a_uz*z/L)**2 - a_vx**2*v_x**2*sin(pi*a_vx*x/L)**2 - a_vz**2*v_z**2*cos(pi*a_vz*z/L)**2 - a_wx**2*w_x**2*cos(pi*a_wx*x/L)**2 -
a_wy**2*w_y**2*cos(pi*a_wy*y/L)**2 - 4*a_ux**2*u_x**2*cos(pi*a_ux*x/L)**2/3 - 4*a_vy**2*v_y**2*cos(pi*a_vy*y/L)**2/3 -
4*a_wz**2*w_z**2*sin(pi*a_wz*z/L)**2/3)/L**2 + pi*DMu_Dx*U*(-4*a_ux*u_x*cos(pi*a_ux*x/L)/3 - 2*a_wz*w_z*sin(pi*a_wz*z/L)/3 +
2*a_vy*v_y*cos(pi*a_vy*y/L)/3)/L + pi*DMu_Dx*V*(a_uy*u_y*sin(pi*a_uy*y/L) + a_vx*v_x*sin(pi*a_vx*x/L))/L + pi*DMu_Dx*W*(a_uz*u_z*sin(pi*a_uz*z/L) -
a_wx*w_x*cos(pi*a_wx*x/L))/L + pi*DMu_Dy*U*(a_uy*u_y*sin(pi*a_uy*y/L) + a_vx*v_x*sin(pi*a_vx*x/L))/L + pi*DMu_Dy*V*(-4*a_vy*v_y*cos(pi*a_vy*y/L)/3 -
2*a_wz*w_z*sin(pi*a_wz*z/L)/3 + 2*a_ux*u_x*cos(pi*a_ux*x/L)/3)/L + pi*DMu_Dy*W*(-a_vz*v_z*cos(pi*a_vz*z/L) - a_wy*w_y*cos(pi*a_wy*y/L))/L +
pi*DMu_Dz*U*(a_uz*u_z*sin(pi*a_uz*z/L) - a_wx*w_x*cos(pi*a_wx*x/L))/L + pi*DMu_Dz*V*(-a_vz*v_z*cos(pi*a_vz*z/L) - a_wy*w_y*cos(pi*a_wy*y/L))/L +
pi*DMu_Dz*W*(2*a_ux*u_x*cos(pi*a_ux*x/L)/3 + 2*a_vy*v_y*cos(pi*a_vy*y/L)/3 + 4*a_wz*w_z*sin(pi*a_wz*z/L)/3)/L + pi*E_t*Rho*(a_ux*u_x*cos(pi*a_ux*x/L)
+ a_vy*v_y*cos(pi*a_vy*y/L) - a_wz*w_z*sin(pi*a_wz*z/L))/L + pi*E_t*a_rhot*rho_t*cos(pi*a_rhot*t/Lt)/Lt + pi*V*a_py*p_y*cos(pi*a_py*y/L)/L +
pi*a_pt*p_t*sin(pi*a_pt*t/Lt)/(Lt*(1 - gamma)) + pi*a_rhox*rho_x*(E_t*U + P*U/(Rho*(1 - gamma)))*cos(pi*a_rhox*x/L)/L + pi*a_rhoy*rho_y*(-E_t*V -
P*V/(Rho*(1 - gamma)))*sin(pi*a_rhoy*y/L)/L + pi*a_rhoz*rho_z*(E_t*W + P*W/(Rho*(1 - gamma)))*cos(pi*a_rhoz*z/L)/L +
kappa*pi**2*(p_x*a_px**2*cos(pi*a_px*x/L) + p_y*a_py**2*sin(pi*a_py*y/L) + p_z*a_pz**2*cos(pi*a_pz*z/L))/(L**2*R*Rho) +
kappa*pi**2*(-2*a_px*a_rhox*p_x*rho_x*cos(pi*a_rhox*x/L)*sin(pi*a_px*x/L) - 2*a_py*a_rhoy*p_y*rho_y*cos(pi*a_py*y/L)*sin(pi*a_rhoy*y/L) -
2*a_pz*a_rhoz*p_z*rho_z*cos(pi*a_rhoz*z/L)*sin(pi*a_pz*z/L))/(L**2*R*Rho**2) - pi*U*a_px*p_x*sin(pi*a_px*x/L)/L - pi*W*a_pz*p_z*sin(pi*a_pz*z/L)/L +
P*kappa*pi**2*(-2*a_rhox**2*rho_x**2*cos(pi*a_rhox*x/L)**2 - 2*a_rhoy**2*rho_y**2*sin(pi*a_rhoy*y/L)**2 -
2*a_rhoz**2*rho_z**2*cos(pi*a_rhoz*z/L)**2)/(L**2*R*Rho**3) + P*kappa*pi**2*(-rho_x*a_rhox**2*sin(pi*a_rhox*x/L) - rho_y*a_rhoy**2*cos(pi*a_rhoy*y/L)
- rho_z*a_rhoz**2*sin(pi*a_rhoz*z/L))/(L**2*R*Rho**2) + pi*P*a_rhot*rho_t*cos(pi*a_rhot*t/Lt)/(Lt*Rho*(1 - gamma)) +
pi*DMu_Dy*a_py*gamma*p_y*cos(pi*a_py*y/L)/(L*Pr*Rho*(1 - gamma)) - pi*DMu_Dx*a_px*gamma*p_x*sin(pi*a_px*x/L)/(L*Pr*Rho*(1 - gamma)) -
pi*DMu_Dz*a_pz*gamma*p_z*sin(pi*a_pz*z/L)/(L*Pr*Rho*(1 - gamma)) + pi*DMu_Dy*P*a_rhoy*gamma*rho_y*sin(pi*a_rhoy*y/L)/(L*Pr*Rho**2*(1 - gamma)) -
pi*DMu_Dx*P*a_rhox*gamma*rho_x*cos(pi*a_rhox*x/L)/(L*Pr*Rho**2*(1 - gamma)) - pi*DMu_Dz*P*a_rhoz*gamma*rho_z*cos(pi*a_rhoz*z/L)/(L*Pr*Rho**2*(1 -
gamma)))
#Q_et= Q_et.subs(kappa,k);
Q_et_num=Q_et.subs(parameters).evalf()  
  
  
print 'Q_et_num-Qet_num=',Q_et_num-Qet_num 
