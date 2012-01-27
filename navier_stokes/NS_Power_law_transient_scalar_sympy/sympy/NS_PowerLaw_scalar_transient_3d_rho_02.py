
# Loading the file for the equations, manufactured solutions and original soource terms

execfile("NS_PowerLaw_scalar_transient_3d_rho_01.py")

# Hierarchic MMS ------------------------------------------------------
var('U,V,W,P,Rho')


# Q1 ------------------------------------------------------------------
Q1n=Q1.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t
* sin(a_rhot * pi * t / Lt),Rho);
Q1n=Q1n.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut *
pi * t / Lt),U);
Q1n=Q1n.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt *
pi * t / Lt),V);
Q1n=Q1n.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt *
pi * t / Lt),W);
Q1n=Q1n.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt *
pi * t / Lt),P);


Q1n=collect(Q1n,Rho*pi/L)


Q_rho_convection=Q1n

# Q2 ------------------------------------------------------------------
Q2n=Q2

Q_rho_time=Q2n


#--------------------------------------------------------------------------
# Factorized source tem for the continuity equation
Q_rho_new=Q1n+Q2n


#--------------------------------------------------------------------------
# Checking factorization

Q=Q_rho_new.subs(Rho,rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t *
sin(a_rhot * pi * t / Lt));
Q=Q.subs(U, u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt));
Q=Q.subs(V, v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt));
Q=Q.subs(W, w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt));
Q=Q.subs(P, p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt));


Res1=expand(Qrho-Q)
Res1==0 #true


Results=Res1
#--------------------------------------------------------------------------
# Writing to C codes
# Saving to C files ---------------------------------------------------------------------
# Asking for futher numerical tests------------------------------------------------------
var('y,Y')
name = input("Save results to C files? [y/n]:") 
if name == y or name==Y:
  if Results == 0:
    #execfile("NS_Power Law_scalar_transient_3d_viscosity.py") #calculating the derivatives of mu
    print 'Successfull factorization: Q_rho'
    print 'Writing C code for Q_rho  -> ../C_codes/NS_PowerLaw_scalar_transient_3d_rho.c'
    execfile("NS_PowerLaw_scalar_transient_3d_rho_codes.py")  
    print 'Done.'
  else:
    print 'ERROR: Possible problems in the factorization!'
    
print '\n'

numerical_tests  = input("Perform numerical test to address correctness in the mapinulations? [y/n]:") 

if numerical_tests == y or numerical_tests == Y:
 execfile("NS_PowerLaw_scalar_transient_3d_rho_check_numerically.py")
 print 'Done.'

