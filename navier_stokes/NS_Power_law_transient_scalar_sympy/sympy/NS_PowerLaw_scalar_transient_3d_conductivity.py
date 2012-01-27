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


#Auxliary calculations for the derivatives of the conductivity kappa

kappa = gamma*R*mu/(gamma-1)/Pr;


# Writing Q_g into a latex file ------------------------------------------



latexQ1=latex(diff(kappa,x))
latexQ2=latex(diff(kappa,y))
latexQ3=latex(diff(kappa,z))

s1=str(latexQ1)
s2=str(latexQ2)
s3=str(latexQ3)

f=open('../latex/derivatives.tex','a')



f.write('\n\n Dkappa_Dx: \n')
f.write(s1)

f.write('\n\n Dkappa_Dy: \n')
f.write(s2)

f.write('\n\n Dkappa_Dz: \n')
f.write(s3)

f.close()


