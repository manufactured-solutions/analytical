#!/usr/bin/env python
# Compute manufactured solution and manufactured forcing at a particular point.
# Used as a test case to ensure C++ implementation consistent with Python, and
# to ensure floating point precision losses are acceptable.

# soln.py contains the nececssary sympy imports
execfile("soln.py")

# Work in (possibly slow) extended precision
PREC = 30;

# Test solution parameters.  These are not physically
# realizable but can be used to check the implementation.
gamma    = 2
R        = 3
beta     = 5
mu_r     = 7
T_r      = 11
k_r      = 13
lambda_r = 17
a_rho0   = 19
a_rhox   = 23
a_rhoxy  = 29
a_rhoxz  = 31
a_rhoy   = 37
a_rhoyz  = 41
a_rhoz   = 43
b_rhox   = 47
b_rhoxy  = 53
b_rhoxz  = 59
b_rhoy   = 61
b_rhoyz  = 67
b_rhoz   = 71
c_rhox   = 73
c_rhoxy  = 79
c_rhoxz  = 83
c_rhoy   = 89
c_rhoyz  = 97
c_rhoz   = 101
d_rhoxy  = 103
d_rhoxz  = 107
d_rhoyz  = 109
e_rhoxy  = 113
e_rhoxz  = 127
e_rhoyz  = 131
f_rho0   = 137
f_rhox   = 139
f_rhoxy  = 149
f_rhoxz  = 151
f_rhoy   = 157
f_rhoyz  = 163
f_rhoz   = 167
g_rho0   = 173
g_rhox   = 179
g_rhoxy  = 181
g_rhoxz  = 191
g_rhoy   = 193
g_rhoyz  = 197
g_rhoz   = 199
a_u0     = 211
a_ux     = 223
a_uxy    = 227
a_uxz    = 229
a_uy     = 233
a_uyz    = 239
a_uz     = 241
b_ux     = 251
b_uxy    = 257
b_uxz    = 263
b_uy     = 269
b_uyz    = 271
b_uz     = 277
c_ux     = 281
c_uxy    = 283
c_uxz    = 293
c_uy     = 307
c_uyz    = 311
c_uz     = 313
d_uxy    = 317
d_uxz    = 331
d_uyz    = 337
e_uxy    = 347
e_uxz    = 349
e_uyz    = 353
f_u0     = 359
f_ux     = 367
f_uxy    = 373
f_uxz    = 379
f_uy     = 383
f_uyz    = 389
f_uz     = 397
g_u0     = 401
g_ux     = 409
g_uxy    = 419
g_uxz    = 421
g_uy     = 431
g_uyz    = 433
g_uz     = 439
a_v0     = 443
a_vx     = 449
a_vxy    = 457
a_vxz    = 461
a_vy     = 463
a_vyz    = 467
a_vz     = 479
b_vx     = 487
b_vxy    = 491
b_vxz    = 499
b_vy     = 503
b_vyz    = 509
b_vz     = 521
c_vx     = 523
c_vxy    = 541
c_vxz    = 547
c_vy     = 557
c_vyz    = 563
c_vz     = 569
d_vxy    = 571
d_vxz    = 577
d_vyz    = 587
e_vxy    = 593
e_vxz    = 599
e_vyz    = 601
f_v0     = 607
f_vx     = 613
f_vxy    = 617
f_vxz    = 619
f_vy     = 631
f_vyz    = 641
f_vz     = 643
g_v0     = 647
g_vx     = 653
g_vxy    = 659
g_vxz    = 661
g_vy     = 673
g_vyz    = 677
g_vz     = 683
a_w0     = 691
a_wx     = 701
a_wxy    = 709
a_wxz    = 719
a_wy     = 727
a_wyz    = 733
a_wz     = 739
b_wx     = 743
b_wxy    = 751
b_wxz    = 757
b_wy     = 761
b_wyz    = 769
b_wz     = 773
c_wx     = 787
c_wxy    = 797
c_wxz    = 809
c_wy     = 811
c_wyz    = 821
c_wz     = 823
d_wxy    = 827
d_wxz    = 829
d_wyz    = 839
e_wxy    = 853
e_wxz    = 857
e_wyz    = 859
f_w0     = 863
f_wx     = 877
f_wxy    = 881
f_wxz    = 883
f_wy     = 887
f_wyz    = 907
f_wz     = 911
g_w0     = 919
g_wx     = 929
g_wxy    = 937
g_wxz    = 941
g_wy     = 947
g_wyz    = 953
g_wz     = 967
a_T0     = 971
a_Tx     = 977
a_Txy    = 983
a_Txz    = 991
a_Ty     = 997
a_Tyz    = 1009
a_Tz     = 1013
b_Tx     = 1019
b_Txy    = 1021
b_Txz    = 1031
b_Ty     = 1033
b_Tyz    = 1039
b_Tz     = 1049
c_Tx     = 1051
c_Txy    = 1061
c_Txz    = 1063
c_Ty     = 1069
c_Tyz    = 1087
c_Tz     = 1091
d_Txy    = 1093
d_Txz    = 1097
d_Tyz    = 1103
e_Txy    = 1109
e_Txz    = 1117
e_Tyz    = 1123
f_T0     = 1129
f_Tx     = 1151
f_Txy    = 1153
f_Txz    = 1163
f_Ty     = 1171
f_Tyz    = 1181
f_Tz     = 1187
g_T0     = 1193
g_Tx     = 1201
g_Txy    = 1213
g_Txz    = 1217
g_Ty     = 1223
g_Tyz    = 1229
g_Tz     = 1231

# Prepare substitution dictionaries for each analytic solution.  Necessary
# since our analytic solution form is repeated for rho, u, v, w, and T.
params_rho = {
    'a_0'  : a_rho0,
    'a_x'  : a_rhox,
    'a_xy' : a_rhoxy,
    'a_xz' : a_rhoxz,
    'a_y'  : a_rhoy,
    'a_yz' : a_rhoyz,
    'a_z'  : a_rhoz,
    'b_x'  : b_rhox,
    'b_xy' : b_rhoxy,
    'b_xz' : b_rhoxz,
    'b_y'  : b_rhoy,
    'b_yz' : b_rhoyz,
    'b_z'  : b_rhoz,
    'c_x'  : c_rhox,
    'c_xy' : c_rhoxy,
    'c_xz' : c_rhoxz,
    'c_y'  : c_rhoy,
    'c_yz' : c_rhoyz,
    'c_z'  : c_rhoz,
    'd_xy' : d_rhoxy,
    'd_xz' : d_rhoxz,
    'd_yz' : d_rhoyz,
    'e_xy' : e_rhoxy,
    'e_xz' : e_rhoxz,
    'e_yz' : e_rhoyz,
    'f_0'  : f_rho0,
    'f_x'  : f_rhox,
    'f_xy' : f_rhoxy,
    'f_xz' : f_rhoxz,
    'f_y'  : f_rhoy,
    'f_yz' : f_rhoyz,
    'f_z'  : f_rhoz,
    'g_0'  : g_rho0,
    'g_x'  : g_rhox,
    'g_xy' : g_rhoxy,
    'g_xz' : g_rhoxz,
    'g_y'  : g_rhoy,
    'g_yz' : g_rhoyz,
    'g_z'  : g_rhoz,
}
params_u = {
    'a_0'  : a_u0,
    'a_x'  : a_ux,
    'a_xy' : a_uxy,
    'a_xz' : a_uxz,
    'a_y'  : a_uy,
    'a_yz' : a_uyz,
    'a_z'  : a_uz,
    'b_x'  : b_ux,
    'b_xy' : b_uxy,
    'b_xz' : b_uxz,
    'b_y'  : b_uy,
    'b_yz' : b_uyz,
    'b_z'  : b_uz,
    'c_x'  : c_ux,
    'c_xy' : c_uxy,
    'c_xz' : c_uxz,
    'c_y'  : c_uy,
    'c_yz' : c_uyz,
    'c_z'  : c_uz,
    'd_xy' : d_uxy,
    'd_xz' : d_uxz,
    'd_yz' : d_uyz,
    'e_xy' : e_uxy,
    'e_xz' : e_uxz,
    'e_yz' : e_uyz,
    'f_0'  : f_u0,
    'f_x'  : f_ux,
    'f_xy' : f_uxy,
    'f_xz' : f_uxz,
    'f_y'  : f_uy,
    'f_yz' : f_uyz,
    'f_z'  : f_uz,
    'g_0'  : g_u0,
    'g_x'  : g_ux,
    'g_xy' : g_uxy,
    'g_xz' : g_uxz,
    'g_y'  : g_uy,
    'g_yz' : g_uyz,
    'g_z'  : g_uz,
}
params_v = {
    'a_0'  : a_v0,
    'a_x'  : a_vx,
    'a_xy' : a_vxy,
    'a_xz' : a_vxz,
    'a_y'  : a_vy,
    'a_yz' : a_vyz,
    'a_z'  : a_vz,
    'b_x'  : b_vx,
    'b_xy' : b_vxy,
    'b_xz' : b_vxz,
    'b_y'  : b_vy,
    'b_yz' : b_vyz,
    'b_z'  : b_vz,
    'c_x'  : c_vx,
    'c_xy' : c_vxy,
    'c_xz' : c_vxz,
    'c_y'  : c_vy,
    'c_yz' : c_vyz,
    'c_z'  : c_vz,
    'd_xy' : d_vxy,
    'd_xz' : d_vxz,
    'd_yz' : d_vyz,
    'e_xy' : e_vxy,
    'e_xz' : e_vxz,
    'e_yz' : e_vyz,
    'f_0'  : f_v0,
    'f_x'  : f_vx,
    'f_xy' : f_vxy,
    'f_xz' : f_vxz,
    'f_y'  : f_vy,
    'f_yz' : f_vyz,
    'f_z'  : f_vz,
    'g_0'  : g_v0,
    'g_x'  : g_vx,
    'g_xy' : g_vxy,
    'g_xz' : g_vxz,
    'g_y'  : g_vy,
    'g_yz' : g_vyz,
    'g_z'  : g_vz,
}
params_w = {
    'a_0'  : a_w0,
    'a_x'  : a_wx,
    'a_xy' : a_wxy,
    'a_xz' : a_wxz,
    'a_y'  : a_wy,
    'a_yz' : a_wyz,
    'a_z'  : a_wz,
    'b_x'  : b_wx,
    'b_xy' : b_wxy,
    'b_xz' : b_wxz,
    'b_y'  : b_wy,
    'b_yz' : b_wyz,
    'b_z'  : b_wz,
    'c_x'  : c_wx,
    'c_xy' : c_wxy,
    'c_xz' : c_wxz,
    'c_y'  : c_wy,
    'c_yz' : c_wyz,
    'c_z'  : c_wz,
    'd_xy' : d_wxy,
    'd_xz' : d_wxz,
    'd_yz' : d_wyz,
    'e_xy' : e_wxy,
    'e_xz' : e_wxz,
    'e_yz' : e_wyz,
    'f_0'  : f_w0,
    'f_x'  : f_wx,
    'f_xy' : f_wxy,
    'f_xz' : f_wxz,
    'f_y'  : f_wy,
    'f_yz' : f_wyz,
    'f_z'  : f_wz,
    'g_0'  : g_w0,
    'g_x'  : g_wx,
    'g_xy' : g_wxy,
    'g_xz' : g_wxz,
    'g_y'  : g_wy,
    'g_yz' : g_wyz,
    'g_z'  : g_wz,
}
params_T = {
    'a_0'  : a_T0,
    'a_x'  : a_Tx,
    'a_xy' : a_Txy,
    'a_xz' : a_Txz,
    'a_y'  : a_Ty,
    'a_yz' : a_Tyz,
    'a_z'  : a_Tz,
    'b_x'  : b_Tx,
    'b_xy' : b_Txy,
    'b_xz' : b_Txz,
    'b_y'  : b_Ty,
    'b_yz' : b_Tyz,
    'b_z'  : b_Tz,
    'c_x'  : c_Tx,
    'c_xy' : c_Txy,
    'c_xz' : c_Txz,
    'c_y'  : c_Ty,
    'c_yz' : c_Tyz,
    'c_z'  : c_Tz,
    'd_xy' : d_Txy,
    'd_xz' : d_Txz,
    'd_yz' : d_Tyz,
    'e_xy' : e_Txy,
    'e_xz' : e_Txz,
    'e_yz' : e_Tyz,
    'f_0'  : f_T0,
    'f_x'  : f_Tx,
    'f_xy' : f_Txy,
    'f_xz' : f_Txz,
    'f_y'  : f_Ty,
    'f_yz' : f_Tyz,
    'f_z'  : f_Tz,
    'g_0'  : g_T0,
    'g_x'  : g_Tx,
    'g_xy' : g_Txy,
    'g_xz' : g_Txz,
    'g_y'  : g_Ty,
    'g_yz' : g_Tyz,
    'g_z'  : g_Tz,
}

# Choose where and when we compute the solution
# This must match the choices made in main.cpp
params_xyzt = {
    'x': Rational(5, 10),
    'y': Rational(6, 10),
    'z': Rational(7, 10),
    't': Rational(8, 10)
}
params_rho.update(params_xyzt)
params_u.update(params_xyzt)
params_v.update(params_xyzt)
params_w.update(params_xyzt)
params_T.update(params_xyzt)

# Compute the analytic solution and derivatives for rho, u, v, w, and T
# These are the necessary inputs to forcing.py which we will source next
rho    = phi   .subs(params_rho)
rho_t  = phi_t .subs(params_rho)
rho_x  = phi_x .subs(params_rho)
rho_xx = phi_xx.subs(params_rho)
rho_xy = phi_xy.subs(params_rho)
rho_xz = phi_xz.subs(params_rho)
rho_y  = phi_y .subs(params_rho)
rho_yy = phi_yy.subs(params_rho)
rho_yz = phi_yz.subs(params_rho)
rho_z  = phi_z .subs(params_rho)
rho_zz = phi_zz.subs(params_rho)

u    = phi   .subs(params_u)
u_t  = phi_t .subs(params_u)
u_x  = phi_x .subs(params_u)
u_xx = phi_xx.subs(params_u)
u_xy = phi_xy.subs(params_u)
u_xz = phi_xz.subs(params_u)
u_y  = phi_y .subs(params_u)
u_yy = phi_yy.subs(params_u)
u_yz = phi_yz.subs(params_u)
u_z  = phi_z .subs(params_u)
u_zz = phi_zz.subs(params_u)

v    = phi   .subs(params_v)
v_t  = phi_t .subs(params_v)
v_x  = phi_x .subs(params_v)
v_xx = phi_xx.subs(params_v)
v_xy = phi_xy.subs(params_v)
v_xz = phi_xz.subs(params_v)
v_y  = phi_y .subs(params_v)
v_yy = phi_yy.subs(params_v)
v_yz = phi_yz.subs(params_v)
v_z  = phi_z .subs(params_v)
v_zz = phi_zz.subs(params_v)

w    = phi   .subs(params_w)
w_t  = phi_t .subs(params_w)
w_x  = phi_x .subs(params_w)
w_xx = phi_xx.subs(params_w)
w_xy = phi_xy.subs(params_w)
w_xz = phi_xz.subs(params_w)
w_y  = phi_y .subs(params_w)
w_yy = phi_yy.subs(params_w)
w_yz = phi_yz.subs(params_w)
w_z  = phi_z .subs(params_w)
w_zz = phi_zz.subs(params_w)

T    = phi   .subs(params_T)
T_t  = phi_t .subs(params_T)
T_x  = phi_x .subs(params_T)
T_xx = phi_xx.subs(params_T)
T_xy = phi_xy.subs(params_T)
T_xz = phi_xz.subs(params_T)
T_y  = phi_y .subs(params_T)
T_yy = phi_yy.subs(params_T)
T_yz = phi_yz.subs(params_T)
T_z  = phi_z .subs(params_T)
T_zz = phi_zz.subs(params_T)

# Compute the forcing using the same equations appearing in writeup.tex
# This also includes term-by-term intermediate values.
execfile("forcing.py")

# Prepare a list of the quantities that we should output.
# This makes assumptions about the names of temporaries in forcing.py,
# but that should not be a big problem.
qoi = """
    rho rho_t rho_x rho_xx rho_xy rho_xz rho_y rho_yy rho_yz rho_z rho_zz
    u   u_t   u_x   u_xx   u_xy   u_xz   u_y   u_yy   u_yz   u_z   u_zz
    v   v_t   v_x   v_xx   v_xy   v_xz   v_y   v_yy   v_yz   v_z   v_zz
    w   w_t   w_x   w_xx   w_xy   w_xz   w_y   w_yy   w_yz   w_z   w_zz
    T   T_t   T_x   T_xx   T_xy   T_xz   T_y   T_yy   T_yz   T_z   T_zz

    gamma  R  beta  mu_r  T_r  k_r  lambda_r

    e e_x e_y e_z e_t
    p p_x p_y p_z mu
    mu_x mu_y mu_z
    lambda_ lambda_x lambda_y lambda_z
    qx_x qy_y qz_z

    rhou rhov rhow
    rhou_x rhov_y rhow_z rhou_t
    rhov_t rhow_t rhoe_t 

    rhouu_x rhouv_y rhouw_z rhouv_x
    rhovv_y rhovw_z rhouw_x rhovw_y
    rhoww_z rhoue_x rhove_y rhowe_z

    tauxx tauyy tauzz tauxy tauxz tauyz
    tauxx_x tauyy_y tauzz_z tauxy_x tauxy_y tauxz_x tauxz_z tauyz_y tauyz_z

    pu_x pv_y pw_z
    utauxx_x vtauxy_x wtauxz_x
    utauxy_y vtauyy_y wtauyz_y
    utauxz_z vtauyz_z wtauzz_z

    Q_rho Q_rhou Q_rhov Q_rhow Q_rhoe
""".split()

# Output C++ namespace containing expected results in as long double constants
def dump(name, value):
    print "    const long double", name, "=",
    if isinstance(value, int):
        print value,
        sys.stdout.write(".0")
    else:
        print value.evalf(PREC),
    sys.stdout.write("L;\n");

print "namespace nsctpl {"
print "namespace test {"
dump('x', params_xyzt['x'])
dump('y', params_xyzt['y'])
dump('z', params_xyzt['z'])
dump('t', params_xyzt['t'])
for quantity in qoi:
    dump(quantity, globals()[quantity])
    sys.stdout.flush()
print "} // end namespace test"
print "} // end namespace nsctpl"
