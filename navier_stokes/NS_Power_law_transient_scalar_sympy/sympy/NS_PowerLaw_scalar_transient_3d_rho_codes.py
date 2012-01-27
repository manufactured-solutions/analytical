
# Writing Q_rho into a C code ------------------------------------------


Q_rho_time=Q2n

Q_rho_convection=Q1n

#unassigning variables Q_rho_time and Q_rho_convection in order to write a more readable C code
var('Q_rho_time')
var('Q_rho_convection')
Q_rho=Q_rho_time+Q_rho_convection

Rho = rho_0 + rho_x * sin(a_rhox * pi * x / L) + rho_y * cos(a_rhoy * pi * y / L) + rho_z * sin(a_rhoz * pi * z / L) + rho_t * sin(a_rhot * pi * t / L);
U = u_0 + u_x * sin(a_ux * pi * x / L) + u_y * cos(a_uy * pi * y / L) + u_z * cos(a_uz * pi * z / L) + u_t * cos(a_ut * pi * t / L);
V = v_0 + v_x * cos(a_vx * pi * x / L) + v_y * sin(a_vy * pi * y / L) + v_z * sin(a_vz * pi * z / L) + v_t * sin(a_vt * pi * t / L);
W = w_0 + w_x * sin(a_wx * pi * x / L) + w_y * sin(a_wy * pi * y / L) + w_z * cos(a_wz * pi * z / L) + w_t * cos(a_wt * pi * t / L);
P = p_0 + p_x * cos(a_px * pi * x / L) + p_y * sin(a_py * pi * y / L) + p_z * cos(a_pz * pi * z / L) + p_t * cos(a_pt * pi * t / L);

from sympy.utilities.codegen import codegen,Routine
codegen((
#	("Rho", Rho),
	("U", U),
	("V", V),
	("W", W),
	("P", P),
	("Q_rho_time",    Q2n   ),
	("Q_rho_convection",    Q1n   ),
        ("Q_rho",    Q_rho  ),

        ), "C", "../C_codes/NS_PowerLaw_scalar_transient_3d_rho", header=True, to_files=True)





# Writing Q_g into a latex file ------------------------------------------
latexQ=latex(Q_rho)


latexQ1=latex(Q1n)#(Q_rho_convection)
latexQ2=latex(Q2n)#(Q_rho_time)

s=str(latexQ)
s1=str(latexQ1)
s2=str(latexQ2)

f=open('../latex/Q_rho.tex','w')

f.write('\n Q_rho: \n')
f.write(s)

f.write('\n\n Q_rho_convection: \n')
f.write(s1)


f.write('\n\n Q_rho_time: \n')
f.write(s2)

#f.write('\n \n')
#f.write(s5)

f.close()




## Writing Q_rho into a file in order to count characters ------------------------------------------
#s=str(Q_rho)
#f=open('../C_code_sympy/SourceQrho_after_factorization.dat','w')
#f.write(s)
#f.close()
