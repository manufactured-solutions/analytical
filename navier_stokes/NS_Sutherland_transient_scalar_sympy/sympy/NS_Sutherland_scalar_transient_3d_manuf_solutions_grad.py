from sympy import *

var('x y z t')

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
phi_an = phi_0 + phi_x * cos(a_phix * pi * x / L) + phi_y * cos(a_phiy * pi * y / L) + phi_z * sin(a_phiz * pi * z / L) + phi_t * cos(a_phit * pi *
t / Lt);


grad_u_an=[diff(u_an,x),diff(u_an,y),diff(u_an,z)]
grad_p_an=[diff(p_an,x),diff(p_an,y),diff(p_an,z)]
grad_phi_an=[diff(phi_an,x),diff(phi_an,y),diff(phi_an,z)]

print(latex(grad_phi_an))

# Writing u_an,v_an, and w_an  and their gradients in a C file ------------------------------------------
dudx=diff(u_an,x);
dudy=diff(u_an,y);
dudz=diff(u_an,z);


dvdx=diff(v_an,x);
dvdy=diff(v_an,y);
dvdz=diff(v_an,z);


dwdx=diff(w_an,x);
dwdy=diff(w_an,y);
dwdz=diff(w_an,z);

dphidx=diff(phi_an,x);
dphidy=diff(phi_an,y);
dphidz=diff(phi_an,z);

drhodx=diff(rho_an,x);
drhody=diff(rho_an,y);
drhodz=diff(rho_an,z);

dpdx=diff(p_an,x);
dpdy=diff(p_an,y);
dpdz=diff(p_an,z);



gradu=ccode(grad_u_an)


from sympy.utilities.codegen import codegen,Routine
codegen((
	
	("u_an", u_an),
	("v_an", v_an),
	("w_an", w_an),
	("p_an", p_an),
	("rho_an", rho_an),
	("phi_an", phi_an),
	("du_dx",dudx),
	("du_dy",dudy),
	("du_dz",dudz),
	("dv_dx",dvdx),
	("dv_dy",dvdy),
	("dv_dz",dvdz),
        ("dw_dx",dwdx),
	("dw_dy",dwdy),
	("dw_dz",dwdz),
	("dp_dx",dpdx),
	("dp_dy",dpdy),
	("dp_dz",dpdz),
	("drho_dx",drhodx),
	("drho_dy",drhody),
	("drho_dz",drhodz),
	("dphi_dx",dphidx),
	("dphi_dy",dphidy),
	("dphi_dz",dphidz),
	
       # ("dw_dz",gradu),
        ), "C", "../C_codes/NS_Sutherland_scalar_transient_manuf_solutions_gradients", header=True, to_files=True)
        



# Writing Q_g into a latex file ------------------------------------------
f=open('../latex/ManufSol_Derivatives.tex','w')
f.write(latex(dudx));
f.write(latex(dudy));
f.write(latex(dudz));
f.write(latex(dvdx));
f.write(latex(dvdy));
f.write(latex(dvdz));
f.write(latex(dwdx));
f.write(latex(dwdy));
f.write(latex(dwdz));
f.write(latex(dpdx));
f.write(latex(dpdy));
f.write(latex(dpdz));
f.write(latex(drhodx));
f.write(latex(drhody));
f.write(latex(drhodz));
f.write(latex(dphidx));
f.write(latex(dphidy));
f.write(latex(dphidz));
f.close()
