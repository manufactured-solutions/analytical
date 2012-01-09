
# Writing Q_w into a C code ------------------------------------------

#Q_w_convection=Q1n
#Q_w_gradp=Q2n
#Q_w_viscous=Q3n
#Q_w_time=Q4n

#unassigning variables in order to write a more readable C code
var('Q_w_time,Q_w_convection,Q_w_gradp,Q_w_viscous')

# Q_w -----------------------------------------------------------------------------------
Q_w=Q_w_convection+Q_w_gradp+Q_w_viscous+Q_w_time

var('T_ref,mu_ref')

Mu = A_mu*T**(3/2)/(T+B_mu);

A_mu = (T_ref + B_mu)* mu_ref/ (T_ref** (3/2)); 

T=P/(R*Rho);

#Rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi*t/Lt);
#U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / Lt);
#V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / Lt);
#W = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / Lt);
#P = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / Lt);

from sympy.utilities.codegen import codegen,Routine
codegen((
	("Rho", rho_an),
	("U", u_an),
	("V", v_an),
	("W", w_an),
	("P", p_an),
	("T", T),
	("A_mu", A_mu),
	("Mu", Mu),
	("DMu_Dx",Dmu_Dx),
	("DMu_Dy",Dmu_Dy),
	("DMu_Dz",Dmu_Dz),
	("Q_w_convection", Q1n),
	("Q_w_gradp", Q2n),
	("Q_w_viscous", Q3n),
	("Q_w_time", Q4n),
	("Q_w",    Q_w   ),

        ), "C", "../C_codes/NS_Sutherland_scalar_transient_3d_w", header=True, to_files=True)


# Writing Q_g into a latex file ------------------------------------------
latexQ=latex(Q_w)
latexQ1=latex(Q1n)#(Q_w_convection)
latexQ2=latex(Q2n)#(Q_w_gradp)
latexQ3=latex(Q3n)#(Q_w_viscous)
latexQ4=latex(Q4n)#(Q_w_time)

s=str(latexQ)
s1=str(latexQ1)
s2=str(latexQ2)
s3=str(latexQ3)
s4=str(latexQ4)
#s5=str('Q_w,Q_w_convection),Q_w_gradp,Q_w_viscous,Q_w_time')

f=open('../latex/Q_w.tex','w')

f.write('\n Q_w: \n')
f.write(s)

f.write('\n\n Q_w_convection: \n')
f.write(s1)

f.write('\n\n Q_w_gradp: \n')
f.write(s2)

f.write('\n\n Q_w_viscous: \n')
f.write(s3)

f.write('\n\n Q_w_time: \n')
f.write(s4)

#f.write('\n \n')
#f.write(s5)

f.close()



## Writing Q_w into a file in order to count characters ------------------------------------------
#s=str(Q_w)
#f=open('../C_code_sympy/SourceQu_after_factorization.dat','w')
#f.write(s)
#f.close()
