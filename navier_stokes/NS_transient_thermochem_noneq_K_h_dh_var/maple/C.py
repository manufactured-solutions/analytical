# This program calculates the source term Q for the 1D NavierStokes equations in themal nonequilibrium (2-temperature model) with chemical reactions N_2 <-> 2N. 
# Two species N and N2
# N -> s1
# N2 -> s2

# STEADY, TWO-TEMPERATURES: T and Tv
# Equilibrium constant K is a to-be-defined function of T: K=K(T).

# First species N density----------------------------------------------------------------------
# rho_N


from sympy import *

var('x t w_dot_N,c_N,w_dot_N2,')

# Auxiliary relations -------------------------------------------------------------------------------------------
#var('R,R_N,R_N2,M_N,M_N2,K,Cf1_N,Cf1_N2,Ea_N,Ea_N2,etaf1_N,etaf1_N2,kf1_N,kf1_N2,kb1_N,kb1_N2')
var('R,R1,M_N,M_N2,Le,Gamma')

var("""L,Lt,rho_N_0,rho_N_x,a_rho_N_x,rho_N_t,a_rho_N_t,rho_N2_0,rho_N2_x,a_rho_N2_x,rho_N2_t,a_rho_N2_t,u_0,u_x,a_ux,u_t,a_ut,
T_0,T_x,a_Tx,T_t,a_Tt,Tv_0,Tv_x,a_Tvx,Tv_t,a_Tvt,""",real=True)

var('RHO_N,RHO_N2,RHO,T,TV,U,M_N,PI')

var('e_elec_N_den,e_elec_N2_den,Sum_eN_thetae2_g_div_e,Sum_eN2_thetae2_g_div_e,E_elec_N,E_elec_N2,E_vib_N2,theta_v_N2,Sum_eN2_thetae3_g_div_e,Sum_eN_thetae3_g_div_e')

var('Mu_N,Mu_N2,Mtot,Phi_N,Phi_N2')

var('A_N,A_N2,B_N,B_N2')

var('Ds,w_dot_V')

var('d_hN2_dT,d_hN_dT,h0_N,h0_N2,Mu_mix,h_N,h_N2')

var('eV_N,eV_N2,eV,DeV_Dx,Kappa_ev_mix,DKappa_ev_Dx,Kappa_tr_mix,DKappa_tr_Dx,Kappa_mix,DKappa_mix_Dx,DCp_Dx,Cp,Cv,DeV_N_Dx,DeV_N2_Dx,DMu_mix_Dx,')
# eV--------------------------------------------

var('Q_eV , Q_eV_time,Q_eV_convection,Q_eV_production,Q_eV_heatflux,Q_eV_diffusion ')

Q_eV = (  Q_eV_time+Q_eV_convection+Q_eV_production+Q_eV_heatflux+Q_eV_diffusion );

Q_eV_time = (  -E_elec_N*a_rho_N_t*PI*rho_N_t*sin(a_rho_N_t*PI*t/Lt)/Lt+(E_vib_N2+E_elec_N2)*a_rho_N2_t*PI*rho_N2_t*cos(a_rho_N2_t*PI*t/Lt)/Lt+a_Tvt*PI*Tv_t*R*RHO_N*cos(a_Tvt*PI*t/Lt)*
Sum_eN_thetae2_g_div_e/Lt/M_N/TV**2/e_elec_N_den+Rational(1,2)*a_Tvt*PI*Tv_t*R*RHO_N2*cos(a_Tvt*PI*t/Lt)*Sum_eN2_thetae2_g_div_e/Lt/M_N/TV**2/e_elec_N2_den+(2*exp(theta_v_N2/TV)*E_vib_N2**2*RHO_N2-RHO_N*E_elec_N**2-2*RHO_N2*E_elec_N2**2)*a_Tvt*PI*Tv_t*M_N*cos(a_Tvt*PI*t/Lt)/Lt/R/TV**2     );


Q_eV_convection = (  eV*a_ux*PI*u_x*RHO*cos(a_ux*PI*x/L)/L+DeV_Dx*RHO*U+(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*eV*PI*U/L     );

Q_eV_production = (  -w_dot_V     );

Q_eV_heatflux = (  a_Tvx**2*PI**2*Tv_x*Kappa_ev_mix*cos(a_Tvx*PI*x/L)/L**2+DKappa_ev_Dx*a_Tvx*PI*Tv_x*sin(a_Tvx*PI*x/L)/L     );

Q_eV_diffusion = (  eV_N*a_rho_N_x**2*PI**2*rho_N_x*Ds*sin(a_rho_N_x*PI*x/L)/L**2+eV_N2*a_rho_N2_x**2*PI**2*rho_N2_x*Ds*cos(a_rho_N2_x*PI*x/L)/L**2+eV_N*DCp_Dx*a_rho_N_x*PI*rho_N_x*
Ds*cos(a_rho_N_x*PI*x/L)/L/Cp-eV_N2*DCp_Dx*a_rho_N2_x*PI*rho_N2_x*Ds*sin(a_rho_N2_x*PI*x/L)/L/Cp-DeV_N_Dx*a_rho_N_x*PI*rho_N_x*Ds*cos(a_rho_N_x*PI*x/L)/L+DeV_N2_Dx*a_rho_N2_x*PI*rho_N2_x*Ds*sin(a_rho_N2_x*PI*x/L)/L+(2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-2*a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*eV_N*a_rho_N_x*PI**2*rho_N_x*Ds*cos(a_rho_N_x*PI*x/L)/L**2/RHO-(2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-2*a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*eV_N2*a_rho_N2_x*PI**2*rho_N2_x*Ds*sin(a_rho_N2_x*PI*x/L)/L**2/RHO-DKappa_mix_Dx*eV_N*a_rho_N_x*PI*rho_N_x*Le*cos(a_rho_N_x*PI*x/L)/L/Cp/RHO+eV_N2*DKappa_mix_Dx*a_rho_N2_x*PI*rho_N2_x*Le*sin(a_rho_N2_x*PI*x/L)/L/Cp/RHO-(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*eV_N*DCp_Dx*PI*Ds*RHO_N/L/Cp/RHO-(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*eV_N2*DCp_Dx*PI*Ds*RHO_N2/L/Cp/RHO+(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*DeV_N_Dx*PI*Ds*RHO_N/L/RHO+(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*DeV_N2_Dx*PI*Ds*RHO_N2/L/RHO-(eV_N*RHO_N+eV_N2*RHO_N2)*(a_rho_N_x**2*rho_N_x*sin(a_rho_N_x*PI*x/L)+a_rho_N2_x**2*rho_N2_x*cos(a_rho_N2_x*PI*x/L))*PI**2*Ds/L**2/RHO-(2*eV_N*RHO_N+2*eV_N2*RHO_N2)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))**2*PI**2*Ds/L**2/RHO**2+(eV_N*RHO_N+eV_N2*RHO_N2)*DKappa_mix_Dx*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*Le/L/Cp/RHO**2     );




 
# E--------------------------------------------
var('Q_E ,  Q_E_time,Q_E_convection,Q_E_work,Q_E_heatflux,Q_E_diffusion    ')

Q_E = (   Q_E_time+Q_E_convection+Q_E_work+Q_E_heatflux+Q_E_diffusion    );
  
Q_E_time = (  
-Rational(3,2)*a_rho_N_t*PI*rho_N_t*R*T*sin(a_rho_N_t*PI*t/Lt)/Lt/M_N+Rational(5,4)*a_rho_N2_t*PI*rho_N2_t*R*T*cos(a_rho_N2_t*PI*t/Lt)/Lt/M_N-a_ut*PI*u_t*RHO*U*sin(a_ut*
PI*t/Lt)/Lt-Rational(1,4)*(6*RHO_N+5*RHO_N2)*a_Tt*PI*T_t*R*sin(a_Tt*PI*t/Lt)/Lt/M_N+a_Tvt*PI*Tv_t*R*RHO_N*cos(a_Tvt*PI*t/Lt)*Sum_eN_thetae2_g_div_e/Lt/M_N/
e_elec_N_den/TV**2+Rational(1,2)*a_Tvt*PI*Tv_t*R*RHO_N2*cos(a_Tvt*PI*t/Lt)*Sum_eN2_thetae2_g_div_e/Lt/M_N/e_elec_N2_den/TV**2-(h0_N+E_elec_N)*a_rho_N_t*PI*
rho_N_t*sin(a_rho_N_t*PI*t/Lt)/Lt+(E_elec_N2+E_vib_N2+h0_N2)*a_rho_N2_t*PI*rho_N2_t*cos(a_rho_N2_t*PI*t/Lt)/Lt-Rational(1,2)*(a_rho_N_t*rho_N_t*sin(a_rho_N_t*PI
*t/Lt)-a_rho_N2_t*rho_N2_t*cos(a_rho_N2_t*PI*t/Lt))*PI*U**2/Lt+2*E_vib_N2**2*M_N*exp(theta_v_N2/TV)*a_Tvt*PI*Tv_t*RHO_N2*cos(a_Tvt*PI*t/Lt)/Lt/R/TV**2
-(RHO_N*E_elec_N**2+2*RHO_N2*E_elec_N2**2)*M_N*a_Tvt*PI*Tv_t*cos(a_Tvt*PI*t/Lt)/Lt/R/TV**2     );

Q_E_convection = (   -2*exp(theta_v_N2/TV)*E_vib_N2**2*a_Tvx*PI*Tv_x*M_N*RHO_N2*U*sin(a_Tvx*PI*x/L)/L/R/TV**2+Rational(5,2)*a_rho_N_x*PI*rho_N_x*R*T*U*cos(a_rho_N_x*PI*x/L)/L/M_N-Rational(7,4)*a_rho_N2_x*PI*rho_N2_x*R*T*U*sin(a_rho_N2_x*PI*x/L)/L/M_N+(h0_N+E_elec_N)*a_rho_N_x*PI*rho_N_x*U*cos(a_rho_N_x*PI*x/L)/L-(E_elec_N2+E_vib_N2+h0_N2)*a_rho_N2_x*PI*rho_N2_x*U*sin(a_rho_N2_x*PI*x/L)/L-Rational(1,4)*(10*RHO_N+7*RHO_N2)*a_Tx*PI*T_x*R*U*sin(a_Tx*PI*x/L)/L/M_N+Rational(1,2)*(2*RHO_N+RHO_N2)*a_ux*PI*u_x*R*T*cos(a_ux*PI*x/L)/L/M_N+(U**2+E)*a_ux*PI*u_x*RHO*cos(a_ux*PI*x/L)/L-a_Tvx*PI*Tv_x*R*RHO_N*U*sin(a_Tvx*PI*x/L)*Sum_eN_thetae2_g_div_e/L/M_N/TV**2/e_elec_N_den-Rational(1,2)*a_Tvx*PI*Tv_x*R*RHO_N2*U*sin(a_Tvx*PI*x/L)*Sum_eN2_thetae2_g_div_e/L/M_N/TV**2/e_elec_N2_den+Rational(1,2)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*U**3/L+(RHO_N*E_elec_N**2+2*RHO_N2*E_elec_N2**2)*a_Tvx*PI*Tv_x*M_N*U*sin(a_Tvx*PI*x/L)/L/R/TV**2     );


Q_E_heatflux = (   a_Tx**2*PI**2*T_x*Kappa_tr_mix*cos(a_Tx*PI*x/L)/L**2+a_Tvx**2*PI**2*Tv_x*Kappa_ev_mix*cos(a_Tvx*PI*x/L)/L**2+DKappa_tr_Dx*a_Tx*PI*T_x*sin(a_Tx*PI*x/L)/L+DKappa_ev_Dx*a_Tvx*PI*Tv_x*sin(a_Tvx*PI*x/L)/L     );

Q_E_work = (  
-Rational(4,3)*a_ux**2*PI**2*u_x**2*Mu_mix*cos(a_ux*PI*x/L)**2/L**2+Rational(4,3)*a_ux**2*PI**2*u_x*Mu_mix*U*sin(a_ux*PI*x/L)/L**2-Rational(4,3)*
DMu_mix_Dx*a_ux*PI*u_x*U*cos(a_ux* PI*x/L)/L     );


Q_E_diffusion = (  
d_hN_dT*a_rho_N_x*a_Tx*PI**2*rho_N_x*T_x*Ds*cos(a_rho_N_x*PI*x/L)*sin(a_Tx*PI*x/L)/L**2-d_hN2_dT*a_rho_N2_x*a_Tx*PI**2*rho_N2_x*T_x*Ds*sin(a_rho_N2_x*
PI*x/L)*sin(a_Tx*PI*x/L)/L**2+h_N*DCp_Dx*a_rho_N_x*PI*rho_N_x*Ds*cos(a_rho_N_x*PI*x/L)/L/Cp-h_N2*DCp_Dx*a_rho_N2_x*PI*rho_N2_x*Ds*sin(a_rho_N2_x*PI*x/
L)/L/Cp+(-a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*d_hN_dT*a_Tx*PI**2*T_x*Ds*RHO_N*sin(a_Tx*PI*x/L)/L**2/
RHO+(-a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*d_hN2_dT*a_Tx*PI**2*T_x*Ds*RHO_N2*sin(a_Tx*PI*x/L)/L**2/RHO-
h_N*DKappa_mix_Dx*a_rho_N_x*PI*rho_N_x*Le*cos(a_rho_N_x*PI*x/L)/L/Cp/RHO-(-2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)+2*a_rho_N2_x*rho_N2_x*sin(
a_rho_N2_x*PI*x/L))*h_N*a_rho_N_x*PI**2*rho_N_x*Ds*cos(a_rho_N_x*PI*x/L)/L**2/RHO+h_N2*DKappa_mix_Dx*a_rho_N2_x*PI*rho_N2_x*Le*sin(a_rho_N2_x*PI*x/L)/
L/Cp/RHO+(-2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)+2*a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*h_N2*a_rho_N2_x*PI**2*rho_N2_x*Ds*sin(a_rho_N2_x*PI
*x/L)/L**2/RHO+(h_N*a_rho_N_x**2*rho_N_x*sin(a_rho_N_x*PI*x/L)+h_N2*a_rho_N2_x**2*rho_N2_x*cos(a_rho_N2_x*PI*x/L))*PI**2*Ds/L**2+(h_N*RHO_N+h_N2*
RHO_N2)*(-a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*DCp_Dx*PI*Ds/L/Cp/RHO-(h_N*RHO_N+h_N2*RHO_N2)*(a_rho_N_x
**2*rho_N_x*sin(a_rho_N_x*PI*x/L)+a_rho_N2_x**2*rho_N2_x*cos(a_rho_N2_x*PI*x/L))*PI**2*Ds/L**2/RHO-(h_N*RHO_N+h_N2*RHO_N2)*(-a_rho_N_x*rho_N_x*cos(
a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*DKappa_mix_Dx*PI*Le/L/Cp/RHO**2-(2*h_N*RHO_N+2*h_N2*RHO_N2)*(-a_rho_N_x*rho_N_x*cos(
a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))**2*PI**2*Ds/L**2/RHO**2     );

  


 
# u--------------------------------------------
var('Q_u , Q_u_time,Q_u_convection,Q_u_gradp,Q_u_work')

Q_u = Q_u_time+Q_u_convection+Q_u_gradp+Q_u_work

Q_u_time = -a_ut*PI*u_t*RHO*sin(a_ut*PI*t/Lt)/Lt-(a_rho_N_t*rho_N_t*sin(a_rho_N_t*PI*t/Lt)-a_rho_N2_t*rho_N2_t*cos(a_rho_N2_t*PI*t/Lt))*PI*U/Lt;

Q_u_convection =2*a_ux*PI*u_x*RHO*U*cos(a_ux*PI*x/L)/L+(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*U**2/L;


Q_u_gradp=-Rational(1,2)*(2*RHO_N+RHO_N2)*a_Tx*PI*T_x*R*sin(a_Tx*PI*x/L)/L/M_N+Rational(1,2)*(2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x
*PI*x/L)) * PI*R*T/L/M_N;


Q_u_work = Rational(4,3)*a_ux**2*PI**2*u_x*Mu_mix*sin(a_ux*PI*x/L)/L**2-Rational(4,3)*DMu_mix_Dx*a_ux*PI*u_x*cos(a_ux*PI*x/L)/L;





# RHO_N--------------------------------------------
var('Q_rho_N,Q_rho_N_time,Q_rho_N_convection,Q_rho_N_diffusion,Q_rho_N_production')
Q_rho_N =Q_rho_N_time+Q_rho_N_convection+Q_rho_N_diffusion+Q_rho_N_production ;

Q_rho_N_time =-a_rho_N_t*PI*rho_N_t*sin(a_rho_N_t*PI*t/Lt)/Lt;

Q_rho_N_convection =a_rho_N_x*PI*rho_N_x*U*cos(a_rho_N_x*PI*x/L)/L+a_ux*PI*u_x*RHO_N*cos(a_ux*PI*x/L)/L;

Q_rho_N_production =-w_dot_N;


Q_rho_N_diffusion =a_rho_N_x**2*PI**2*rho_N_x*Ds*sin(a_rho_N_x*PI*x/L)/L**2+DCp_Dx*a_rho_N_x*PI*rho_N_x*Le*Kappa_mix*cos(a_rho_N_x*PI*x/L)/L/Cp**2/RHO-(-a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*a_rho_N_x*PI**2*rho_N_x*Ds*cos(a_rho_N_x*PI*x/L)/L**2/RHO-(-a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*a_rho_N_x*PI**2*rho_N_x*Le*Kappa_mix*cos(a_rho_N_x*PI*x/L)/L**2/Cp/RHO**2-DKappa_mix_Dx*a_rho_N_x*PI*rho_N_x*Le*cos(a_rho_N_x*PI*x/L)/L/Cp/RHO-(a_rho_N_x**2*rho_N_x*sin(a_rho_N_x*PI*x/L)+a_rho_N2_x**2*rho_N2_x*cos(a_rho_N2_x*PI*x/L))*PI**2*Ds*RHO_N/L**2/RHO+(-a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*DCp_Dx*PI*Le*Kappa_mix*RHO_N/L/Cp**2/RHO**2-(-a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))**2*PI**2*Ds*RHO_N/L**2/RHO**2-(-a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))**2*PI**2*Le*Kappa_mix*RHO_N/L**2/Cp/RHO**3-(-a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*DKappa_mix_Dx*PI*Le*RHO_N/L/Cp/RHO**2;




# RHO_N2--------------------------------------------
var('Q_rho_N2,Q_rho_N2_time,Q_rho_N2_convection,Q_rho_N2_diffusion,Q_rho_N2_production')

Q_rho_N2 =Q_rho_N2_time+Q_rho_N2_convection+Q_rho_N2_diffusion+Q_rho_N2_production;

Q_rho_N2_time =a_rho_N2_t*PI*rho_N2_t*cos(a_rho_N2_t*PI*t/Lt)/Lt;

Q_rho_N2_convection =-a_rho_N2_x*PI*rho_N2_x*U*sin(a_rho_N2_x*PI*x/L)/L+a_ux*PI*u_x*RHO_N2*cos(a_ux*PI*x/L)/L;

Q_rho_N2_production =-w_dot_N2;

Q_rho_N2_diffusion=(
a_rho_N2_x**2*PI**2*rho_N2_x*Ds*cos(a_rho_N2_x*PI*x/L)/L**2-DCp_Dx*a_rho_N2_x*PI*rho_N2_x*Le*Kappa_mix*sin(a_rho_N2_x*PI*x/L)/L/Cp**2/RHO+(-a_rho_N_x
*rho_N_x*cos(a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*a_rho_N2_x*PI**2*rho_N2_x*Ds*sin(a_rho_N2_x*PI*x/L)/L**2/RHO+(-a_rho_N_x*
rho_N_x*cos(a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*a_rho_N2_x*PI**2*rho_N2_x*Le*Kappa_mix*sin(a_rho_N2_x*PI*x/L)/L**2/Cp/RHO**2
+DKappa_mix_Dx*a_rho_N2_x*PI*rho_N2_x*Le*sin(a_rho_N2_x*PI*x/L)/L/Cp/RHO-(a_rho_N_x**2*rho_N_x*sin(a_rho_N_x*PI*x/L)+a_rho_N2_x**2*rho_N2_x*cos(
a_rho_N2_x*PI*x/L))*PI**2*Ds*RHO_N2/L**2/RHO+(-a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*DCp_Dx*PI*Le*
Kappa_mix*RHO_N2/L/Cp**2/RHO**2-(-a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))**2*PI**2*Ds*RHO_N2/L**2/RHO**2-(
-a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))**2*PI**2*Le*Kappa_mix*RHO_N2/L**2/Cp/RHO**3-(-a_rho_N_x*rho_N_x*
cos(a_rho_N_x*PI*x/L)+a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*DKappa_mix_Dx*PI*Le*RHO_N2/L/Cp/RHO**2  );





#=====================================================


Cp = (  Cv+R)   ; 

Cv=sympify(Rational(3,2)*RHO_N*R/RHO/M_N+Rational(5,4)*RHO_N2*R/RHO/M_N+RHO_N*R*Sum_eN_thetae2_g_div_e/RHO/TV**2/M_N/e_elec_N_den-RHO_N*E_elec_N**2*M_N/RHO/R/TV**2+2*
RHO_N2*exp(theta_v_N2/TV)*M_N*E_vib_N2**2/RHO/R/TV**2+Rational(1,2)*RHO_N2*R*Sum_eN2_thetae2_g_div_e/RHO/TV**2/M_N/e_elec_N2_den-2*RHO_N2*E_elec_N2**2*M_N/RHO/R
/TV**2 ) 



Kappa_mix = (   Kappa_ev_mix+Kappa_tr_mix     );


Kappa_ev_mix = (  Mu_N*R*RHO_N*Mtot*Sum_eN_thetae2_g_div_e/
M_N**2/RHO/Phi_N/TV**2/e_elec_N_den-Mu_N*RHO_N*Mtot*E_elec_N**2/RHO/Phi_N/R/TV**2+
Rational(1,4)*Mu_N2*R*RHO_N2*Mtot*Sum_eN2_thetae2_g_div_e/M_N**2/RHO/Phi_N2/TV**2/
e_elec_N2_den-Mu_N2*RHO_N2*Mtot*E_elec_N2**2/RHO/Phi_N2/R/TV**2+Mu_N2*exp(
theta_v_N2/TV)*RHO_N2*E_vib_N2**2*Mtot/RHO/Phi_N2/R/TV**2     );


Kappa_tr_mix = (   Rational(15,4)*Mu_N*RHO_N*Mtot*R/RHO/M_N**2/Phi_N+Rational(19,16)*Mu_N2*RHO_N2*Mtot*R/RHO/M_N**2/Phi_N2     );


DMu_mix_Dx= (   -(2*A_N*ln(T)-2*A_N2*ln(T)+B_N-B_N2)*a_Tx*PI*T_x*Mu_N*Mtot*RHO_N*sin(a_Tx*PI*x/L)/L/M_N/Phi_N/RHO/T-Rational(1,24)*sqrt(Rational(6,1))*(1+Rational(1,2)*sqrt(Mu_N2/Mu_N)*Rational(2,1)**Rational(3,4))**2*a_rho_N_x*PI*rho_N_x*Mu_N2*Mtot**2*RHO_N2*cos(a_rho_N_x*PI*x/L)/L/M_N**2/Phi_N2**2/RHO**2+Rational(1,12)*sqrt(Rational(3,1))*(1+sqrt(Mu_N/Mu_N2)*Rational(2,1)**Rational(1,4))**2*a_rho_N2_x*PI*rho_N2_x*Mu_N*Mtot**2*RHO_N*sin(a_rho_N2_x*PI*x/L)/L/M_N**2/Phi_N**2/RHO**2+Rational(1,12)*sqrt(Rational(3,1))*Rational(2,1)**Rational(1,4)*(1+sqrt(Mu_N/Mu_N2)*Rational(2,1)**Rational(1,4))*(2*A_N*ln(T)-2*A_N2*ln(T)+B_N-B_N2)*a_Tx*PI*T_x*Mu_N**2*Mtot**2*RHO_N*RHO_N2*sin(a_Tx*PI*x/L)/sqrt(Mu_N/Mu_N2)/L/M_N**2/Mu_N2/Phi_N**2/RHO**2/T-Rational(1,48)*sqrt(Rational(6,1))*Rational(2,1)**Rational(3,4)*(1+Rational(1,2)*sqrt(Mu_N2/Mu_N)*Rational(2,1)**Rational(3,4))*(2*A_N*ln(T)-2*A_N2*ln(T)+B_N-B_N2)*a_Tx*PI*T_x*Mu_N2**2*Mtot**2*RHO_N*RHO_N2*sin(a_Tx*PI*x/L)/sqrt(Mu_N2/Mu_N)/L/M_N**2/Mu_N/Phi_N2**2/RHO**2/T+(RHO*M_N*Phi_N-RHO_N*Mtot)*a_rho_N_x*PI*rho_N_x*Mu_N*Mtot*cos(a_rho_N_x*PI*x/L)/L/M_N**2/Phi_N**2/RHO**2-Rational(1,4)*(2*RHO*M_N*Phi_N2-RHO_N2*Mtot)*a_rho_N2_x*PI*rho_N2_x*Mu_N2*Mtot*sin(a_rho_N2_x*PI*x/L)/L/M_N**2/Phi_N2**2/RHO**2-Rational(1,2)*(2*Mu_N*Phi_N2*RHO_N+Mu_N2*Phi_N*RHO_N2)*(2*A_N2*ln(T)+B_N2)*a_Tx*PI*T_x*Mtot*sin(a_Tx*PI*x/L)/L/M_N/Phi_N/Phi_N2/RHO/T+Rational(1,24)*sqrt(Rational(3,1))*(2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(1+sqrt(Mu_N/Mu_N2)*Rational(2,1)**Rational(1,4))**2*PI*Mu_N*Mtot**3*RHO_N*RHO_N2/L/M_N**3/Phi_N**2/RHO**3+Rational(1,48)*sqrt(Rational(6,1))*(2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(1+Rational(1,2)*sqrt(Mu_N2/Mu_N)*Rational(2,1)**Rational(3,4))**2*PI*Mu_N2*Mtot**3*RHO_N*RHO_N2/L/M_N**3/Phi_N2**2/RHO**3-Rational(1,4)*(2*Mu_N*Phi_N2*RHO_N+Mu_N2*Phi_N*RHO_N2)*(2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*Mtot**2/L/M_N**2/Phi_N/Phi_N2/RHO**2+Rational(1,24)*sqrt(Rational(3,1))*(2*RHO*M_N-2*RHO_N*Mtot-RHO_N2*Mtot)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(1+sqrt(Mu_N/Mu_N2)*Rational(2,1)**Rational(1,4))**2*PI*Mu_N*Mtot**2*RHO_N*RHO_N2/L/M_N**3/Phi_N**2/RHO**4+Rational(1,48)*sqrt(Rational(6,1))*(2*RHO*M_N-2*RHO_N*Mtot-RHO_N2*Mtot)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(1+Rational(1,2)*sqrt(Mu_N2/Mu_N)*Rational(2,1)**Rational(3,4))**2*PI*Mu_N2*Mtot**2*RHO_N*RHO_N2/L/M_N**3/Phi_N2**2/RHO**4+Rational(1,8)*(4*Mu_N*Phi_N2**2*RHO_N**2+Mu_N2*Phi_N**2*RHO_N2**2)*(2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*Mtot**3/L/M_N**3/Phi_N**2/Phi_N2**2/RHO**3-Rational(1,8)*(2*RHO*M_N-2*RHO_N*Mtot-RHO_N2*Mtot)*(4*M_N*Mu_N*Phi_N*Phi_N2**2*RHO_N*RHO+2*M_N*Mu_N2*Phi_N**2*Phi_N2*RHO_N2*RHO-4*Mu_N*Phi_N2**2*Mtot*RHO_N**2-Mu_N2*Phi_N**2*Mtot*RHO_N2**2)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*Mtot/L/M_N**3/Phi_N**2/Phi_N2**2/RHO**4     );

DKappa_mix_Dx = (   DKappa_ev_Dx+ DKappa_tr_Dx     );

DKappa_ev_Dx = (  -(2*Mu_N*Phi_N2*E_elec_N**2*RHO_N+2*Mu_N2*Phi_N*E_elec_N2**2*RHO_N2)*a_Tvx*PI*Tv_x*Mtot*sin(a_Tvx*PI*x/L)/L/R/Phi_N/Phi_N2/RHO/TV**3- 4*exp(theta_v_N2/TV)**2*a_Tvx*PI*Tv_x*M_N*Mu_N2*Mtot*E_vib_N2**3*RHO_N2*sin(a_Tvx*PI*x/L)/L/R**2/Phi_N2/RHO/TV**4-Rational(1,4)*(2*M_N*Phi_N2*RHO-RHO_N2*Mtot)*(2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*exp(theta_v_N2/TV)*PI*Mu_N2*Mtot**2*E_vib_N2**2*RHO_N2/L/M_N**2/R/Phi_N2**2/RHO**3/TV**2 -(M_N*Phi_N*RHO-RHO_N*Mtot)*a_rho_N_x*PI*rho_N_x*Mu_N*Mtot*E_elec_N**2*cos(a_rho_N_x*PI*x/L)/L/M_N/R/Phi_N**2/RHO**2/TV**2+Rational(1,2)*(2*M_N*Phi_N2*RHO-RHO_N2*Mtot)*a_rho_N2_x*PI*rho_N2_x*Mu_N2*Mtot*E_elec_N2**2*sin(a_rho_N2_x*PI*x/L)/L/M_N/R/Phi_N2**2/RHO**2/TV**2+(2*A_N*ln(T)-2*A_N2*ln(T)+B_N-B_N2)*a_Tx*PI*T_x*Mu_N*Mtot*E_elec_N**2*RHO_N*sin(a_Tx*PI*x/L)/L/R/Phi_N/RHO/T/TV**2-(2*A_N2*ln(T)+B_N2)*a_Tx*PI*T_x*R*Mu_N*Mtot*RHO_N*sin(a_Tx*PI*x/L)*Sum_eN_thetae2_g_div_e/L/M_N**2/Phi_N/RHO/T/TV**2/e_elec_N_den-(2*A_N*ln(T)-2*A_N2*ln(T)+B_N-B_N2)*a_Tx*PI*T_x*R*Mu_N*Mtot*RHO_N*sin(a_Tx*PI*x/L)*Sum_eN_thetae2_g_div_e/L/M_N**2/Phi_N/RHO/T/TV**2/e_elec_N_den+Rational(1,12)*sqrt(Rational(3,1))*(1+sqrt(Mu_N/Mu_N2)*Rational(2,1)**Rational(1,4))**2*a_rho_N2_x*PI*rho_N2_x*R*Mu_N*Mtot**2*RHO_N*sin(a_rho_N2_x*PI*x/L)*Sum_eN_thetae2_g_div_e/L/M_N**3/Phi_N**2/RHO**2/TV**2/e_elec_N_den+Rational(1,12)*sqrt(Rational(3,1))*Rational(2,1)**Rational(1,4)*(1+sqrt(Mu_N/Mu_N2)*Rational(2,1)**Rational(1,4))*(2*A_N*ln(T)-2*A_N2*ln(T)+B_N-B_N2)*a_Tx*PI*T_x*R*Mu_N**2*Mtot**2*RHO_N*RHO_N2*sin(a_Tx*PI*x/L)*Sum_eN_thetae2_g_div_e/sqrt(Mu_N/Mu_N2)/L/M_N**3/Mu_N2/Phi_N**2/RHO**2/T/TV**2/e_elec_N_den+(3*M_N*E_elec_N+2*R*TV)*a_Tvx*PI*Tv_x*Mu_N*Mtot*RHO_N*sin(a_Tvx*PI*x/L)*Sum_eN_thetae2_g_div_e/L/M_N**2/Phi_N/RHO/TV**4/e_elec_N_den+(M_N*Phi_N*RHO-RHO_N*Mtot)*a_rho_N_x*PI*rho_N_x*R*Mu_N*Mtot*cos(a_rho_N_x*PI*x/L)*Sum_eN_thetae2_g_div_e/L/M_N**3/Phi_N**2/RHO**2/TV**2/e_elec_N_den+Rational(1,24)*sqrt(Rational(3,1))*(2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(1+sqrt(Mu_N/Mu_N2)*Rational(2,1)**Rational(1,4))**2*PI*R*Mu_N*Mtot**3*RHO_N*RHO_N2*Sum_eN_thetae2_g_div_e/L/M_N**4/Phi_N**2/RHO**3/TV**2/e_elec_N_den+Rational(1,24)*sqrt(Rational(3,1))*(2*RHO*M_N-2*RHO_N*Mtot-RHO_N2*Mtot)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(1+sqrt(Mu_N/Mu_N2)*Rational(2,1)**Rational(1,4))**2*PI*R*Mu_N*Mtot**2*RHO_N*RHO_N2*Sum_eN_thetae2_g_div_e/L/M_N**4/Phi_N**2/RHO**4/TV**2/e_elec_N_den-Rational(1,2)*(M_N*Phi_N*RHO-RHO_N*Mtot)*(2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*R*Mu_N*Mtot**2*RHO_N*Sum_eN_thetae2_g_div_e/L/M_N**4/Phi_N**2/RHO**3/TV**2/e_elec_N_den-Rational(1,2)*(2*RHO*M_N-2*RHO_N*Mtot-RHO_N2*Mtot)*(M_N*Phi_N*RHO-RHO_N*Mtot)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*R*Mu_N*Mtot*RHO_N*Sum_eN_thetae2_g_div_e/L/M_N**4/Phi_N**2/RHO**4/TV**2/e_elec_N_den-Rational(1,4)*Sum_eN2_thetae3_g_div_e*a_Tvx*PI*Tv_x*R*Mu_N2*Mtot*RHO_N2*sin(a_Tvx*PI*x/L)/L/M_N**2/Phi_N2/RHO/TV**4/e_elec_N2_den-Rational(1,4)*(2*A_N2*ln(T)+B_N2)*a_Tx*PI*T_x*R*Mu_N2*Mtot*RHO_N2*sin(a_Tx*PI*x/L)*Sum_eN2_thetae2_g_div_e/L/M_N**2/Phi_N2/RHO/T/TV**2/e_elec_N2_den- Rational(1,48)*sqrt(Rational(6,1))*(1+Rational(1,2)*sqrt(Mu_N2/Mu_N)*Rational(2,1)**Rational(3,4))**2*a_rho_N_x*PI*rho_N_x*R*Mu_N2*Mtot**2*RHO_N2*cos(a_rho_N_x*PI*x/L)*Sum_eN2_thetae2_g_div_e/L/M_N**3/Phi_N2**2/RHO**2/TV**2/e_elec_N2_den-Rational(1,96)*sqrt(Rational(6,1))*Rational(2,1)**Rational(3,4)*(1+Rational(1,2)*sqrt(Mu_N2/Mu_N)*Rational(2,1)**Rational(3,4))*(2*A_N*ln(T)-2*A_N2*ln(T)+B_N-B_N2)*a_Tx*PI*T_x*R*Mu_N2**2*Mtot**2*RHO_N*RHO_N2*sin(a_Tx*PI*x/L)*Sum_eN2_thetae2_g_div_e/sqrt(Mu_N2/Mu_N)/L/M_N**3/Mu_N/Phi_N2**2/RHO**2/T/TV**2/e_elec_N2_den+Rational(1,2)*(3*M_N*E_elec_N2+R*TV)*a_Tvx*PI*Tv_x*Mu_N2*Mtot*RHO_N2*sin(a_Tvx*PI*x/L)*Sum_eN2_thetae2_g_div_e/L/M_N**2/Phi_N2/RHO/TV**4/e_elec_N2_den-Rational(1,8)*(2*M_N*Phi_N2*RHO-RHO_N2*Mtot)*a_rho_N2_x*PI*rho_N2_x*R*Mu_N2*Mtot*sin(a_rho_N2_x*PI*x/L)*Sum_eN2_thetae2_g_div_e/L/M_N**3/Phi_N2**2/RHO**2/TV**2/e_elec_N2_den+(Mu_N*Phi_N2*E_elec_N**2*RHO_N+Mu_N2*Phi_N*E_elec_N2**2*RHO_N2-Mu_N2*exp(theta_v_N2/TV)*Phi_N*E_vib_N2**2*RHO_N2)*(2*A_N2*ln(T)+B_N2)*a_Tx*PI*T_x*Mtot*sin(a_Tx*PI*x/L)/L/R/Phi_N/Phi_N2/RHO/T/TV**2-Rational(1,12)*sqrt(Rational(3,1))*Rational(2,1)**Rational(1,4)*(1+sqrt(Mu_N/Mu_N2)*Rational(2,1)**Rational(1,4))*(2*A_N*ln(T)-2*A_N2*ln(T)+B_N-B_N2)*a_Tx*PI*T_x*Mu_N**2*Mtot**2*E_elec_N**2*RHO_N*RHO_N2*sin(a_Tx*PI*x/L)/sqrt(Mu_N/Mu_N2)/L/M_N/R/Mu_N2/Phi_N**2/RHO**2/T/TV**2+theta_v_N2*exp(theta_v_N2/TV)*a_Tvx*PI*Tv_x*Mu_N2*Mtot*E_vib_N2**2*RHO_N2*sin(a_Tvx*PI*x/L)/L/R/Phi_N2/RHO/TV**4-Rational(1,12)*sqrt(Rational(6,1))*(-E_elec_N2**2+exp(theta_v_N2/TV)*E_vib_N2**2)*(1+Rational(1,2)*sqrt(Mu_N2/Mu_N)*Rational(2,1)**Rational(3,4))**2*a_rho_N_x*PI*rho_N_x*Mu_N2*Mtot**2*RHO_N2*cos(a_rho_N_x*PI*x/L)/L/M_N/R/Phi_N2**2/RHO**2/TV**2-Rational(1,2)*(2*M_N*Phi_N2*RHO-RHO_N2*Mtot)*exp(theta_v_N2/TV)*a_rho_N2_x*PI*rho_N2_x*Mu_N2*Mtot*E_vib_N2**2*sin(a_rho_N2_x*PI*x/L)/L/M_N/R/Phi_N2**2/RHO**2/TV**2-Rational(1,24)*sqrt(Rational(6,1))*Rational(2,1)**Rational(3,4)*(-E_elec_N2**2+exp(theta_v_N2/TV)*E_vib_N2**2)*(1+Rational(1,2)*sqrt(Mu_N2/Mu_N)*Rational(2,1)**Rational(3,4))*(2*A_N*ln(T)-2*A_N2*ln(T)+B_N-B_N2)*a_Tx*PI*T_x*Mu_N2**2*Mtot**2*RHO_N*RHO_N2*sin(a_Tx*PI*x/L)/sqrt(Mu_N2/Mu_N)/L/M_N/R/Mu_N/Phi_N2**2/RHO**2/T/TV**2-(2*Mu_N*Phi_N2*E_elec_N**3*RHO_N+4*Mu_N2*Phi_N*E_elec_N2**3*RHO_N2)*a_Tvx*PI*Tv_x*M_N*Mtot*sin(a_Tvx*PI*x/L)/L/R**2/Phi_N/Phi_N2/RHO/TV**4-Rational(1,24)*sqrt(Rational(3,1))*(2*RHO*M_N-2*RHO_N*Mtot-RHO_N2*Mtot)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(1+sqrt(Mu_N/Mu_N2)*Rational(2,1)**Rational(1,4))**2*PI*Mu_N*Mtot**2*E_elec_N**2*RHO_N*RHO_N2/L/M_N**2/R/
Phi_N**2/RHO**4/TV**2+Rational(1,24)*sqrt(Rational(6,1))*(-E_elec_N2**2+exp(theta_v_N2/TV)*E_vib_N2**2)*(2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(1+Rational(1,2)*sqrt(Mu_N2/Mu_N)*Rational(2,1)**Rational(3,4))**2*PI*Mu_N2*Mtot**3*RHO_N*RHO_N2/L/M_N**2/R/Phi_N2**2/RHO**3/TV**2+Rational(1,24)*sqrt(Rational(6,1))*(2*RHO*M_N-2*RHO_N*Mtot-RHO_N2*Mtot)*(-E_elec_N2**2+exp(theta_v_N2/TV)*E_vib_N2**2)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(1+Rational(1,2)*sqrt(Mu_N2/Mu_N)*Rational(2,1)**Rational(3,4))**2*PI*Mu_N2*Mtot**2*RHO_N*RHO_N2/L/M_N**2/R/Phi_N2**2/RHO**4/TV**2-Rational(1,4)*(2*RHO*M_N-2*RHO_N*Mtot-RHO_N2*Mtot)*(2*M_N*Phi_N2*RHO-RHO_N2*Mtot)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*exp(theta_v_N2/TV)*PI*Mu_N2*Mtot*E_vib_N2**2*RHO_N2/L/M_N**2/R/Phi_N2**2/RHO**4/TV**2-Rational(1,24)*sqrt(Rational(3,1))*(-4*sqrt(Rational(3,1))*M_N*Mu_N*Phi_N*Phi_N2**2*E_elec_N**2*RHO_N*RHO-4*sqrt(Rational(3,1))*M_N*Mu_N2*Phi_N**2*Phi_N2*E_elec_N2**2*RHO_N2*RHO+4*sqrt(Rational(3,1))*Mu_N*Phi_N2**2*Mtot*E_elec_N**2*RHO_N**2+(1+sqrt(Mu_N/Mu_N2)*Rational(2,1)**Rational(1,4))**2*Mu_N*Phi_N2**2*Mtot*E_elec_N**2*RHO_N*RHO_N2+2*sqrt(Rational(3,1))*Mu_N2*Phi_N**2*Mtot*E_elec_N2**2*RHO_N2**2)*(2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*Mtot**2/L/M_N**2/R/Phi_N**2/Phi_N2**2/RHO**3/TV**2+Rational(1,4)*(2*RHO*M_N-2*RHO_N*Mtot-RHO_N2*Mtot)*(2*M_N*Mu_N*Phi_N*Phi_N2**2*E_elec_N**2*RHO_N*RHO+2*M_N*Mu_N2*Phi_N**2*Phi_N2*E_elec_N2**2*RHO_N2*RHO-2*Mu_N*Phi_N2**2*Mtot*E_elec_N**2*RHO_N**2-Mu_N2*Phi_N**2*Mtot*E_elec_N2**2*RHO_N2**2)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*Mtot/L/M_N**2/R/Phi_N**2/Phi_N2**2/RHO**4/TV**2+2*exp(theta_v_N2/TV)*a_Tvx*PI*Tv_x*Mu_N2*Mtot*E_vib_N2**2*RHO_N2*sin(a_Tvx*PI*x/L)/L/R/Phi_N2/RHO/TV**3-Sum_eN_thetae3_g_div_e*a_Tvx*PI*Tv_x*R*Mu_N*Mtot*RHO_N*sin(a_Tvx*PI*x/L)/L/M_N**2/Phi_N/RHO/TV**4/e_elec_N_den+Rational(1,96)*sqrt(Rational(6,1))*(2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(1+Rational(1,2)*sqrt(Mu_N2/Mu_N)*Rational(2,1)**Rational(3,4))**2*PI*R*Mu_N2*Mtot**3*RHO_N*RHO_N2*Sum_eN2_thetae2_g_div_e/L/M_N**4/Phi_N2**2/RHO**3/TV**2/e_elec_N2_den+Rational(1,96)*sqrt(Rational(6,1))*(2*RHO*M_N-2*RHO_N*Mtot-RHO_N2*Mtot)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(1+Rational(1,2)*sqrt(Mu_N2/Mu_N)*Rational(2,1)**Rational(3,4))**2*PI*R*Mu_N2*Mtot**2*RHO_N*RHO_N2*Sum_eN2_thetae2_g_div_e/L/M_N**4/Phi_N2**2/RHO**4/TV**2/e_elec_N2_den-Rational(1,16)*(2*M_N*Phi_N2*RHO-RHO_N2*Mtot)*(2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*R*Mu_N2*Mtot**2*RHO_N2*Sum_eN2_thetae2_g_div_e/L/M_N**4/Phi_N2**2/RHO**3/TV**2/e_elec_N2_den-Rational(1,16)*(2*RHO*M_N-2*RHO_N*Mtot-RHO_N2*Mtot)*(2*M_N*Phi_N2*RHO-RHO_N2*Mtot)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*R*Mu_N2*Mtot*RHO_N2*Sum_eN2_thetae2_g_div_e/L/M_N**4/Phi_N2**2/RHO**4/TV**2/e_elec_N2_den-Rational(1,12)*sqrt(Rational(3,1))*(1+sqrt(Mu_N/Mu_N2)*Rational(2,1)**Rational(1,4))**2*a_rho_N2_x*PI*rho_N2_x*Mu_N*Mtot**2*E_elec_N**2*RHO_N*sin(a_rho_N2_x*PI*x/L)/L/M_N/R/Phi_N**2/RHO**2/TV**2     );

DKappa_tr_Dx = (  
-Rational(15,4)*(2*A_N*ln(T)-2*A_N2*ln(T)+B_N-B_N2)*a_Tx*PI*T_x*R*Mu_N*Mtot*RHO_N*sin(a_Tx*PI*x/L)/L/M_N**2/Phi_N/RHO/T-Rational(19,192)*sqrt(Rational(6,1))*(1+Rational(1,2)*sqrt(Mu_N2/Mu_N)*Rational(2,1)**Rational(3,4))**2*a_rho_N_x*PI*rho_N_x*R*Mu_N2*Mtot**2*RHO_N2*cos(a_rho_N_x*PI*x/L)/L/M_N**3/Phi_N2**2/RHO**2+Rational(5,16)*sqrt(Rational(3,1))*(1+sqrt(Mu_N/Mu_N2)*Rational(2,1)**Rational(1,4))**2*
a_rho_N2_x*PI*rho_N2_x*R*Mu_N*Mtot**2*RHO_N*sin(a_rho_N2_x*PI*x/L)/L/M_N**3/Phi_N**2/RHO**2+Rational(5,16)*sqrt(Rational(3,1))*Rational(2,1)**Rational(1,4)*(1+sqrt(Mu_N/Mu_N2)*Rational(2,1)**Rational(1,4))*(2*A_N
*ln(T)-2*A_N2*ln(T)+B_N-B_N2)*a_Tx*PI*T_x*R*Mu_N**2*Mtot**2*RHO_N*RHO_N2*sin(a_Tx*PI*x/L)/sqrt(Mu_N/Mu_N2)/L/M_N**3/Mu_N2/Phi_N**2/RHO**2/T-Rational(19,384)*
sqrt(Rational(6,1))*Rational(2,1)**Rational(3,4)*(1+Rational(1,2)*sqrt(Mu_N2/Mu_N)*Rational(2,1)**Rational(3,4))*(2*A_N*ln(T)-2*A_N2*ln(T)+B_N-B_N2)*a_Tx*PI*T_x*R*Mu_N2**2*Mtot**2*RHO_N*RHO_N2*sin(a_Tx*PI*x/L)/
sqrt(Mu_N2/Mu_N)/L/M_N**3/Mu_N/Phi_N2**2/RHO**2/T+Rational(15,4)*(RHO*M_N*Phi_N-RHO_N*Mtot)*a_rho_N_x*PI*rho_N_x*R*Mu_N*Mtot*cos(a_rho_N_x*PI*x/L)/L/M_N**3/
Phi_N**2/RHO**2-Rational(19,32)*(2*RHO*M_N*Phi_N2-RHO_N2*Mtot)*a_rho_N2_x*PI*rho_N2_x*R*Mu_N2*Mtot*sin(a_rho_N2_x*PI*x/L)/L/M_N**3/Phi_N2**2/RHO**2-Rational(1,16)*(60*
Mu_N*Phi_N2*RHO_N+19*Mu_N2*Phi_N*RHO_N2)*(2*A_N2*ln(T)+B_N2)*a_Tx*PI*T_x*R*Mtot*sin(a_Tx*PI*x/L)/L/M_N**2/Phi_N/Phi_N2/RHO/T+Rational(5,32)*sqrt(Rational(3,1))*(2*a_rho_N_x
*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(1+sqrt(Mu_N/Mu_N2)*Rational(2,1)**Rational(1,4))**2*PI*R*Mu_N*Mtot**3*RHO_N*RHO_N2/L/M_N**4/
Phi_N**2/RHO**3+Rational(19,384)*sqrt(Rational(6,1))*(2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(1+Rational(1,2)*sqrt(Mu_N2/Mu_N)*Rational(2,1)**Rational(3,4)
)**2*PI*R*Mu_N2*Mtot**3*RHO_N*RHO_N2/L/M_N**4/Phi_N2**2/RHO**3+Rational(5,32)*sqrt(Rational(3,1))*(2*RHO*M_N-2*RHO_N*Mtot-RHO_N2*Mtot)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x
/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(1+sqrt(Mu_N/Mu_N2)*Rational(2,1)**Rational(1,4))**2*PI*R*Mu_N*Mtot**2*RHO_N*RHO_N2/L/M_N**4/Phi_N**2/RHO**4+Rational(19,384)*sqrt(6
)*(2*RHO*M_N-2*RHO_N*Mtot-RHO_N2*Mtot)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)
- a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(1+Rational(1,2)*sqrt(Mu_N2/Mu_N)*Rational(2,1)**Rational(3,4))**2*PI*R*Mu_N2*Mtot**2*RHO_N*RHO_N2/L/M_N**4/Phi_N2**2/RHO**4 -
Rational(15,8)*(2*RHO*M_N-2*RHO_N*Mtot-RHO_N2*Mtot)*(RHO*M_N*Phi_N-RHO_N*Mtot)*(a_rho_N_x*
rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*R*Mu_N*Mtot*RHO_N/L/M_N**4/Phi_N**2/RHO**4-Rational(19,64)*(2*RHO*M_N-2*RHO_N*Mtot-
RHO_N2*Mtot)*(2*RHO*M_N*Phi_N2-RHO_N2*Mtot)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*R*Mu_N2*Mtot*
RHO_N2/L/M_N**4/Phi_N2**2/RHO**4-Rational(1,64)*(120*M_N*Mu_N*Phi_N*Phi_N2**2*RHO_N*RHO+38*M_N*Mu_N2*Phi_N**2*Phi_N2*RHO_N2*RHO-120*Mu_N*Phi_N2**2*Mtot*RHO_N**2
-19*Mu_N2*Phi_N**2*Mtot*RHO_N2**2)*(2*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-
a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*R*Mtot**2/L/M_N**4/Phi_N**2/Phi_N2**2/RHO**3     );



DCp_Dx = ( 
2*a_Tvx*PI*Tv_x*M_N*theta_v_N2*exp(theta_v_N2/TV)*E_vib_N2**2*RHO_N2*sin(a_Tvx*PI*x/L)/L/R/RHO/TV**4-a_rho_N_x*PI*rho_N_x*M_N*E_elec_N**2*cos(
a_rho_N_x*PI*x/L)/L/R/RHO/TV**2+a_rho_N_x*PI*rho_N_x*R*cos(a_rho_N_x*PI*x/L)*Sum_eN_thetae2_g_div_e/L/M_N/RHO/TV**2/e_elec_N_den-(-2*E_elec_N2**2+2*
exp(theta_v_N2/TV)*E_vib_N2**2)*a_rho_N2_x*PI*rho_N2_x*M_N*sin(a_rho_N2_x*PI*x/L)/L/R/RHO/TV**2-Rational(1,2)*a_rho_N2_x*PI*rho_N2_x*R*sin(a_rho_N2_x*PI*x/L)*
Sum_eN2_thetae2_g_div_e/L/M_N/RHO/TV**2/e_elec_N2_den+2*a_Tvx*PI*Tv_x*R*RHO_N*sin(a_Tvx*PI*x/L)*Sum_eN_thetae2_g_div_e/L/M_N/RHO/TV**3/e_elec_N_den+
a_Tvx*PI*Tv_x*R*RHO_N2*sin(a_Tvx*PI*x/L)*Sum_eN2_thetae2_g_div_e/L/M_N/RHO/TV**3/e_elec_N2_den+3*a_Tvx*PI*Tv_x*E_elec_N*RHO_N*sin(a_Tvx*PI*x/L)*
Sum_eN_thetae2_g_div_e/L/RHO/TV**4/e_elec_N_den+3*a_Tvx*PI*Tv_x*E_elec_N2*RHO_N2*sin(a_Tvx*PI*x/L)*Sum_eN2_thetae2_g_div_e/L/RHO/TV**4/e_elec_N2_den+(
4*exp(theta_v_N2/TV)*E_vib_N2**2*RHO_N2-2*RHO_N*E_elec_N**2-4*RHO_N2*E_elec_N2**2)*a_Tvx*PI*Tv_x*M_N*sin(a_Tvx*PI*x/L)/L/R/RHO/TV**3-a_Tvx*PI*Tv_x*R*
RHO_N*sin(a_Tvx*PI*x/L)*Sum_eN_thetae3_g_div_e/L/M_N/RHO/TV**4/e_elec_N_den-Rational(1,2)*a_Tvx*PI*Tv_x*R*RHO_N2*sin(a_Tvx*PI*x/L)*Sum_eN2_thetae3_g_div_e/L/M_N
/RHO/TV**4/e_elec_N2_den+Rational(1,4)*(6*a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-5*a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*R/L/M_N/RHO-(8*exp(theta_v_N2
/TV)**2*E_vib_N2**3*RHO_N2+2*E_elec_N**3*RHO_N+8*E_elec_N2**3*RHO_N2)*a_Tvx*PI*Tv_x*M_N**2*sin(a_Tvx*PI*x/L)/L/R**2/RHO/TV**4-Rational(1,4)*(a_rho_N_x*rho_N_x*
cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(6*RHO_N+5*RHO_N2)*PI*R/L/M_N/RHO**2-(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-
a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*R*RHO_N*Sum_eN_thetae2_g_div_e/L/M_N/RHO**2/TV**2/e_elec_N_den-Rational(1,2)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*
x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*R*RHO_N2*Sum_eN2_thetae2_g_div_e/L/M_N/RHO**2/TV**2/e_elec_N2_den-(a_rho_N_x*rho_N_x*cos(a_rho_N_x
*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*(2*exp(theta_v_N2/TV)*E_vib_N2**2*RHO_N2-RHO_N*E_elec_N**2
-2*RHO_N2*E_elec_N2**2)*PI*M_N/L/R/RHO**2/TV**2     );


DeV_Dx = (  DeV_N_Dx+DeV_N2_Dx     );

DeV_N_Dx = (  a_rho_N_x*PI*rho_N_x*E_elec_N*cos(a_rho_N_x*PI*x/L)/L/RHO+a_Tvx*PI*Tv_x*M_N*E_elec_N**2*RHO_N*sin(a_Tvx*PI*x/L)/L/R/RHO/TV**2-a_Tvx*PI*Tv_x*R*RHO_N*sin(a_Tvx*PI*x/L)*Sum_eN_thetae2_g_div_e/L/M_N/RHO/TV**2/e_elec_N_den-(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*E_elec_N*RHO_N/L/RHO**2     );

  
DeV_N2_Dx = (  -(E_vib_N2+E_elec_N2)*a_rho_N2_x*PI*rho_N2_x*sin(a_rho_N2_x*PI*x/L)/L/RHO-(-2*E_elec_N2**2+2*exp(theta_v_N2/TV)*E_vib_N2**2)*a_Tvx*PI*Tv_x*M_N*RHO_N2*sin(a_Tvx*PI*x/L)/L/R/RHO/TV**2-Rational(1,2)*a_Tvx*PI*Tv_x*R*RHO_N2*sin(a_Tvx*PI*x/L)*Sum_eN2_thetae2_g_div_e/L/M_N/RHO/TV**2/e_elec_N2_den-(E_vib_N2+E_elec_N2)*(a_rho_N_x*rho_N_x*cos(a_rho_N_x*PI*x/L)-a_rho_N2_x*rho_N2_x*sin(a_rho_N2_x*PI*x/L))*PI*RHO_N2/L/RHO**2     );



# Writing source terms into latex form ---------------------------------------
# Rho_N ----------------------------------------------------------------------
latexQ=latex(Q_rho_N)
latexQ1=latex(Q_rho_N_convection)
latexQ2=latex(Q_rho_N_diffusion)
latexQ3=latex(Q_rho_N_production)
latexQ5=latex(Q_rho_N_time)

s=str(latexQ)
s1=str(latexQ1)
s2=str(latexQ2)
s3=str(latexQ3)
s5=str(latexQ5)

f=open('../latex/Source_latex_02.tex','w')

f.write('\n Q_rho_N: \n')
f.write(s)

f.write('\n\n Q_rho_N_convection: \n')
f.write(s1)

f.write('\n\n Q_rho_N_diffusion: \n')
f.write(s2)

f.write('\n\n Q_rho_N_production: \n')
f.write(s3)

f.write('\n\n Q_rho_N_time: \n')
f.write(s5)

f.write('\n---------------------------------------- \n')

# Rho_N2 ----------------------------------------------------------------------
latexQ=latex(Q_rho_N2)
latexQ1=latex(Q_rho_N2_convection)
latexQ2=latex(Q_rho_N2_diffusion)
latexQ3=latex(Q_rho_N2_production)
latexQ5=latex(Q_rho_N2_time)

s=str(latexQ)
s1=str(latexQ1)
s2=str(latexQ2)
s3=str(latexQ3)
s5=str(latexQ5)

f.write('\n Q_rho_N2: \n')
f.write(s)

f.write('\n\n Q_rho_N2_convection: \n')
f.write(s1)

f.write('\n\n Q_rho_N2_diffusion: \n')
f.write(s2)

f.write('\n\n Q_rho_N2_production: \n')
f.write(s3)

f.write('\n\n Q_rho_N2_time: \n')
f.write(s5)

f.write('\n---------------------------------------- \n')

# u ----------------------------------------------------------------------
latexQ=latex(Q_u)
latexQ1=latex(Q_u_convection)
latexQ2=latex(Q_u_gradp)
latexQ3=latex(Q_u_work)
latexQ5=latex(Q_u_time)

s=str(latexQ)
s1=str(latexQ1)
s2=str(latexQ2)
s3=str(latexQ3)
s5=str(latexQ5)

f.write('\n Q_u: \n')
f.write(s)

f.write('\n\n Q_u_convection: \n')
f.write(s1)

f.write('\n\n Q_u_gradp: \n')
f.write(s2)

f.write('\n\n Q_u_work: \n')
f.write(s3)

f.write('\n\n Q_u_time: \n')
f.write(s5)

f.write('\n---------------------------------------- \n')

# E ----------------------------------------------------------------------
latexQ=latex(Q_E)
latexQ1=latex(Q_E_convection)
latexQ2=latex(Q_E_diffusion)
latexQ3=latex(Q_E_work)
latexQ4=latex(Q_E_heatflux)
latexQ5=latex(Q_E_time)

s=str(latexQ)
s1=str(latexQ1)
s2=str(latexQ2)
s3=str(latexQ3)
s4=str(latexQ4)
s5=str(latexQ5)

f.write('\n Q_et: \n')
f.write(s)

f.write('\n\n Q_E_convection: \n')
f.write(s1)

f.write('\n\n Q_E_diffusion: \n')
f.write(s2)

f.write('\n\n Q_E_work: \n')
f.write(s3)

f.write('\n\n Q_E_heatflux: \n')
f.write(s4)

f.write('\n\n Q_E_time: \n')
f.write(s5)

f.write('\n---------------------------------------- \n')

# eV----------------------------------------------------
latexQ=latex(Q_eV)
latexQ1=latex(Q_eV_convection)
latexQ2=latex(Q_eV_diffusion)
latexQ3=latex(Q_eV_production)
latexQ4=latex(Q_eV_heatflux)
latexQ5=latex(Q_eV_time)

s=str(latexQ)
s1=str(latexQ1)
s2=str(latexQ2)
s3=str(latexQ3)
s4=str(latexQ4)
s5=str(latexQ5)


f.write('\n Q_eV: \n')
f.write(s)

f.write('\n\n Q_eV_convection: \n')
f.write(s1)

f.write('\n\n Q_eV_diffusion: \n')
f.write(s2)

f.write('\n\n Q_eV_production: \n')
f.write(s3)

f.write('\n\n Q_eV_heatflux: \n')
f.write(s4)

f.write('\n\n Q_eV_time: \n')
f.write(s5)

f.write('\n---------------------------------------- \n')



f.write('\n Cv: \n')
f.write( latex(Cv) ) 


f.write('\n Kappa_mix : \n')
f.write( latex(Kappa_mix ) ) 

f.write('\n Kappa_ev_mix: \n')
f.write( latex(Kappa_ev_mix) ) 

f.write('\n Kappa_tr_mix: \n')
f.write( latex(Kappa_tr_mix) ) 

f.write('\n DMu_mix_Dx: \n')
f.write( latex(DMu_mix_Dx) ) 

f.write('\n DKappa_mix_Dx : \n')
f.write( latex(DKappa_mix_Dx) ) 

f.write('\n DKappa_tr_Dx : \n')
f.write( latex(DKappa_tr_Dx  ) ) 

f.write('\n DKappa_ev_Dx : \n')
f.write( latex(DKappa_ev_Dx ) ) 

f.write('\n DCp_Dx : \n')
f.write( latex(DCp_Dx ) ) 

f.write('\n DeV_Dx : \n')
f.write( latex(DeV_Dx ) ) 

f.write('\n DeV_N_Dx: \n')
f.write( latex(DeV_N_Dx) ) 


f.write('\n DeV_N2_Dx: \n')
f.write( latex(DeV_N2_Dx) ) 

f.close()





