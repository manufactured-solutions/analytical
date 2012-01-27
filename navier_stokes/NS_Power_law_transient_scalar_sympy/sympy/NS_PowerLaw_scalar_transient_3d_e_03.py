
# Loading the file for the equations, manufactured solutions and original soource terms

execfile("NS_PowerLaw_scalar_transient_3d_e_02.py")



Q_et=Q1n+Q2n+Q3n+Q4n+Q5n
Q=Q_et.subs(E_t,-P/(Rho*(1 - gamma)) + U**2/2 + V**2/2 + W**2/2)
Q=Q.subs(Rho,rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot *
pi * t / Lt));
Q=Q.subs(U, u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt));
Q=Q.subs(V, v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt));
Q=Q.subs(W, w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt));
Q=Q.subs(P, p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt));
Q=Q.subs(Mu,mu);
Q=Q.subs(DMu_Dx,diff(mu,x));
Q=Q.subs(DMu_Dy,diff(mu,y));
Q=Q.subs(DMu_Dz,diff(mu,z));

Q_et=Q

Qet=Q1+Q2+Q3+Q4+Q5;

Res=expand(Q_et-Qet)
Res==0 #false! :(

# The expression is too complex for Sympy and it gets lost
Res=(Q_et-Qet)
Res=Res.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot *
pi * t / Lt),Rho);
Res=Res.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt),U);
Res=Res.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt),V);
Res=Res.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt),W);
Res=Res.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt),P);
Res=expand(Res)

Res==0 #true

var('Results')

if Res == 0 :
  execfile("NS_PowerLaw_scalar_transient_3d_viscosity.py") #calculating the derivatives of mu
  print 'Successfull factorization: Q_et'
  #print 'Writing C code for Q_et  -> ../C_codes/NS_PowerLaw_scalar_transient_3d_e.c'
  #execfile("NS_PowerLaw_scalar_transient_3d_e_codes.py")  
  #print 'Done.'
else:
  # Sympy can't handle it, so let's check the terms separately:
  #
  #Q1 and Q1n -----------------------------------------------------------------------------
  print 'Checking Q1'
  Q=Q1n.subs(E_t,-P/(Rho*(1 - gamma)) + U**2/2 + V**2/2 + W**2/2)
  Q=Q.subs(Rho,rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot *
  pi * t / Lt));
  Q=Q.subs(U, u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt));
  Q=Q.subs(V, v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt));
  Q=Q.subs(W, w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt));
  Q=Q.subs(P, p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt));
  Q=Q.subs(Mu,mu);
  Q=Q.subs(DMu_Dx,diff(mu,x));
  Q=Q.subs(DMu_Dy,diff(mu,y));
  Q=Q.subs(DMu_Dz,diff(mu,z));
  #
  Res1=expand(Q1-Q)
  Res1==0 #false
  #
  # The expression is too complex and Sympy gets lost
  Res1=Q1-Q;
  Res1=Res1.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot *
pi * t / Lt),Rho);
  Res1=Res1.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt),U);
  Res1=Res1.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt),V);
  Res1=Res1.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt),W);
  Res1=Res1.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt),P);
  Res1=expand(Res1)
  Res1==0 #true
  #
  #Q2 and Q2n -----------------------------------------------------------------------------
  print 'Checking Q2'
  Q=Q2n.subs(E_t,-P/(Rho*(1 - gamma)) + U**2/2 + V**2/2 + W**2/2)
  Q=Q.subs(Rho,rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot *
  pi * t / Lt));
  Q=Q.subs(U, u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt));
  Q=Q.subs(V, v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt));
  Q=Q.subs(W, w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt));
  Q=Q.subs(P, p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt));
  Q=Q.subs(Mu,mu);
  Q=Q.subs(DMu_Dx,diff(mu,x));
  Q=Q.subs(DMu_Dy,diff(mu,y));
  Q=Q.subs(DMu_Dz,diff(mu,z));
  #
  Res2=expand(Q2-Q)
  Res2==0 #true
  #
  #Q3 and Q3n -----------------------------------------------------------------------------
  print 'Checking Q3'
  Q=Q3n.subs(E_t,-P/(Rho*(1 - gamma)) + U**2/2 + V**2/2 + W**2/2)
  Q=Q.subs(Rho,rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot *
pi * t / Lt));
  Q=Q.subs(U, u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt));
  Q=Q.subs(V, v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt));
  Q=Q.subs(W, w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt));
  Q=Q.subs(P, p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt));
  Q=Q.subs(LAMBDA, (alpha-Rational(2,3))*Mu);
  Q=Q.subs(Mu,mu);
  Q=Q.subs(DMu_Dx,diff(mu,x));
  Q=Q.subs(DMu_Dy,diff(mu,y));
  Q=Q.subs(DMu_Dz,diff(mu,z));
  #
  Res3=expand(Q3-Q)
  Res3==0 #true
  #
  #Q4 and Q4n -----------------------------------------------------------------------------
  print 'Checking Q4'
  #Q=Q4n.subs(E_t,-P/(Rho*(1 - gamma)) + U**2/2 + V**2/2 + W**2/2)
  #Q=Q.subs(Rho,rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot *
#pi * t / Lt));
  #Q=Q.subs(U, u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt));
  #Q=Q.subs(V, v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt));
  #Q=Q.subs(W, w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt));
  #Q=Q.subs(P, p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt));
  #Q=Q.subs(Mu,mu);
  #Q=Q.subs(DMu_Dx,diff(mu,x));
  #Q=Q.subs(DMu_Dy,diff(mu,y));
  #Q=Q.subs(DMu_Dz,diff(mu,z));
  ##
  #Res4=expand(Q4-Q)
  #Res4==0 #true
  # It doesn't work otherwise
  Res4=Q4-Q4n;
  Res4=Res4.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot *
pi * t / Lt),Rho);
  Res4=Res4.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt),U);
  Res4=Res4.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt),V);
  Res4=Res4.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt),W);
  Res4=Res4.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt),P);
  Res4=Res4.subs(diff(k,x),-R*gamma*DMu_Dx/(1 - gamma)/Pr)
  Res4=Res4.subs(diff(k,y),-R*gamma*DMu_Dy/(1 - gamma)/Pr)
  Res4=Res4.subs(diff(k,z),-R*gamma*DMu_Dz/(1 - gamma)/Pr)
  Res4=Res4.subs(k,kappa)
  Res4=expand(Res4)
  Res4==0 
  #
  #Q5 and Q5n -----------------------------------------------------------------------------
  print 'Checking Q5'
  Q=Q5n.subs(E_t,-P/(Rho*(1 - gamma)) + U**2/2 + V**2/2 + W**2/2)
  Q=Q.subs(Rho,rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot *
pi * t / Lt));
  Q=Q.subs(U, u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt));
  Q=Q.subs(V, v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt));
  Q=Q.subs(W, w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt));
  Q=Q.subs(P, p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt));
  Q=Q.subs(Mu,mu);
  Q=Q.subs(DMu_Dx,diff(mu,x));
  Q=Q.subs(DMu_Dy,diff(mu,y));
  Q=Q.subs(DMu_Dz,diff(mu,z));
  #
  Res5=expand(Q5-Q)
  Res5==0 #true
  #
 # The expression is too complex and Sympy gets lost
  Res5=Q5-Q;
  Res5=Res5.subs(rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot *
pi * t / Lt),Rho);
  Res5=Res5.subs( u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t /
Lt),U);
  Res5=Res5.subs( v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt),V);
  Res5=Res5.subs( w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt),W);
  Res5=Res5.subs( p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt),P);
  Res5=expand(Res5)
  #
  Res5==0 #true
  #
  #------------------------------------------------------------------------------------------ 
  Results=Res1+Res2+Res3+Res4+Res5
  if Results == 0: 
    execfile("NS_PowerLaw_scalar_transient_3d_viscosity.py") #calculating the derivatives of mu
    print 'Successfull factorization: Q_et'
    #print 'Writing C code for Q_et  -> ../C_codes/NS_PowerLaw_scalar_transient_3d_e.c'
    #execfile("NS_PowerLaw_scalar_transient_3d_e_codes.py")  
    #print 'Done.'
  else:
    print 'ERROR: Possible problems in the factorization!'
  

 
  
 
 
 
# Saving to C files ---------------------------------------------------------------------
# Asking for futher numerical tests------------------------------------------------------
var('y,Y')

print '\n'

name = input("Save results to C files? [y/n]:") 
if name == y or name==Y:
  if Results == 0 or Res==0:
    print 'Writing C code for Q_et  -> ../C_codes/NS_PowerLaw_scalar_transient_3d_e.c'
    execfile("NS_PowerLaw_scalar_transient_3d_e_codes.py")  
    print 'Done.'
    
print '\n'

numerical_tests  = input("Perform numerical test to address correctness in the mapinulations? [y/n]:") 

if numerical_tests == y or numerical_tests == Y:
 execfile("NS_PowerLaw_scalar_transient_3d_e_check_numerically.py")
 print 'Done.'

