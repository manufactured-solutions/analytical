
# Writing Q_rho into a C code ------------------------------------------


Q_rho_time=Q2n

Q_rho_convection=Q1n

#unassigning variables Q_rho_time and Q_rho_convection in order to write a more readable C code
var('Q_rho_time')
var('Q_rho_convection')
Q_rho=Q_rho_time+Q_rho_convection


from sympy.utilities.codegen import codegen,Routine
codegen((
	("U", u_an),
	("V", v_an),
	("W", w_an),
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
