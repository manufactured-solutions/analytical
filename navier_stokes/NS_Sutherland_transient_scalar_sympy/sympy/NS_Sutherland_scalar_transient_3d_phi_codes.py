
# Writing Q_phi into a C code ------------------------------------------


#Q_phi_time=Q2n
#Q_phi_convection=Q1n
#Q_phi_diffusion=Q3n

#unassigning variables Q_phi_time and Q_phi_convection in order to write a more readable C code
var('Q_phi_time, Q_phi_convection, Q_phi_diffusion')

Q_phi=Q_phi_time+Q_phi_convection+Q_phi_diffusion


from sympy.utilities.codegen import codegen,Routine
codegen((
#	("Rho", Rho),
	("U", u_an),
	("V", v_an),
	("W", w_an),
	#("P", p_an),
	("Phi", phi_an),
	("Q_phi_time",    Q2n   ),
	("Q_phi_convection",    Q1n   ),
	("Q_phi_diffusion", Q3n ),
        ("Q_phi",    Q_phi  ),

        ), "C", "../C_codes/NS_Sutherland_scalar_transient_3d_phi", header=True, to_files=True)





# Writing Q_g into a latex file ------------------------------------------
latexQ=latex(Q_phi)


latexQ1=latex(Q1n)#(Q_phi_convection)
latexQ2=latex(Q2n)#(Q_phi_time)
latexQ3=latex(Q3n)#(Q_phi_diffusion)

s=str(latexQ)
s1=str(latexQ1)
s2=str(latexQ2)
s3=str(latexQ3)

f=open('../latex/Q_phi.tex','w')

f.write('\n Q_phi: \n')
f.write(s)

f.write('\n\n Q_phi_convection: \n')
f.write(s1)

f.write('\n\n Q_phi_time: \n')
f.write(s2)

f.write('\n\n Q_phi_diffusion: \n')
f.write(s3)

f.close()




## Writing Q_phi into a file in order to count characters ------------------------------------------
#s=str(Q_phi)
#f=open('../C_code_sympy/SourceQrho_after_factorization.dat','w')
#f.write(s)
#f.close()
