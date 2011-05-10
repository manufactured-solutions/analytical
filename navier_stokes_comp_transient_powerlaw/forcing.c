/* Assuming that we are given                                  */
/*     rho, rho_t, rho_x, rho_xy, rho_xz, rho_y, rho_yz, rho_z */
/*     u,   u_t,   u_x,   u_xy,   u_xz,   u_y,   u_yz,   u_z   */
/*     v,   v_t,   v_x,   v_xy,   v_xz,   v_y,   v_yz,   v_z   */
/*     w,   w_t,   w_x,   w_xy,   w_xz,   w_y,   w_yz,   w_z   */
/*     T,   T_t,   T_x,   T_xy,   T_xz,   T_y,   T_yz,   T_z   */
/* and the coefficients                                        */
/*     gamma, R, beta, mu_r, T_r, k_r, lambda_r                */
/* compute the source terms                                    */
/*     Q_rho, Q_rhou, Q_rhov, Q_rhow, Q_rhoe                   */
/* necessary to force the solution rho, u, v, w, and T.        */

/* Computations stemming from the constitutive relationships */
e        = R * T   / (gamma - 1) + (u*u   + v*v   + w*w  ) / 2;
e_x      = R * T_x / (gamma - 1) + (u*u_x + v*v_x + w*w_x);
e_y      = R * T_y / (gamma - 1) + (u*u_y + v*v_y + w*w_y);
e_z      = R * T_z / (gamma - 1) + (u*u_z + v*v_z + w*w_z);
e_t      = R * T_t / (gamma - 1) + (u*u_t + v*v_t + w*w_t);
p        = rho * R * T;
p_x      = rho_x * R * T + rho * R * T_x;
p_y      = rho_y * R * T + rho * R * T_y;
p_z      = rho_z * R * T + rho * R * T_z;
mu       = mu_r * pow(T / T_r, beta);
mu_x     = beta * mu_r / T_r * pow(T / T_r, beta - 1) * T_x;
mu_y     = beta * mu_r / T_r * pow(T / T_r, beta - 1) * T_y;
mu_z     = beta * mu_r / T_r * pow(T / T_r, beta - 1) * T_z;
lambda   = lambda_r / mu_r * mu;
lambda_x = lambda_r / mu_r * mu_x;
lambda_y = lambda_r / mu_r * mu_y;
lambda_z = lambda_r / mu_r * mu_z;
qx_x     = k_r / mu_r * (mu_x * T_x + mu * T_xx);
qy_y     = k_r / mu_r * (mu_y * T_y + mu * T_yy);
qz_z     = k_r / mu_r * (mu_z * T_z + mu * T_zz);

/* Computations stemming from the compressible, Newtonian fluid model */
rhou    = rho * u;
rhov    = rho * v;
rhow    = rho * w;
rhou_x  = rho_x * u + rho * u_x;
rhov_y  = rho_y * v + rho * v_y;
rhow_z  = rho_z * w + rho * w_z;
rhou_t  = rho_t * u + rho * u_t;
rhov_t  = rho_t * v + rho * v_t;
rhow_t  = rho_t * w + rho * w_t;
rhoe_t  = rho_t * e + rho * e_t;

rhouu_x = (rho_x * u * u) + (rho * u_x * u) + (rho * u * u_x);
rhouv_y = (rho_y * u * v) + (rho * u_y * v) + (rho * u * v_y);
rhouw_z = (rho_z * u * w) + (rho * u_z * w) + (rho * u * w_z);
rhouv_x = (rho_x * u * v) + (rho * u_x * v) + (rho * u * v_x);
rhovv_y = (rho_y * v * v) + (rho * v_y * v) + (rho * v * v_y);
rhovw_z = (rho_z * v * w) + (rho * v_z * w) + (rho * v * w_z);
rhouw_x = (rho_x * u * w) + (rho * u_x * w) + (rho * u * w_x);
rhovw_y = (rho_y * v * w) + (rho * v_y * w) + (rho * v * w_y);
rhoww_z = (rho_z * w * w) + (rho * w_z * w) + (rho * w * w_z);
rhoue_x = (rho_x * u * e) + (rho * u_x * e) + (rho * u * e_x);
rhove_y = (rho_y * v * e) + (rho * v_y * e) + (rho * v * e_y);
rhowe_z = (rho_z * w * e) + (rho * w_z * e) + (rho * w * e_z);

tauxx = mu * (u_x + u_x) + lambda * (u_x + v_y + w_z);
tauyy = mu * (v_y + v_y) + lambda * (u_x + v_y + w_z);
tauzz = mu * (w_z + w_z) + lambda * (u_x + v_y + w_z);
tauxy = mu * (u_y + v_x);
tauxz = mu * (u_z + w_x);
tauyz = mu * (v_z + w_y);

tauxx_x = mu_x * (u_x  + u_x ) + lambda_x * (u_x  + v_y  + w_z )
        + mu   * (u_xx + u_xx) + lambda   * (u_xx + v_xy + w_xz);
tauyy_y = mu_y * (v_y  + v_y ) + lambda_y * (u_x  + v_y  + w_z )
        + mu   * (v_yy + v_yy) + lambda   * (u_xy + v_yy + w_yz);
tauzz_z = mu_z * (w_z  + w_z ) + lambda_z * (u_x  + v_y  + w_z )
        + mu   * (w_zz + w_zz) + lambda   * (u_xz + v_yz + w_zz);

tauxy_x = mu_x * (u_y + v_x) + mu * (u_xy + v_xx);
tauxy_y = mu_y * (u_y + v_x) + mu * (u_yy + v_xy);
tauxz_x = mu_x * (u_z + w_x) + mu * (u_xz + w_xx);
tauxz_z = mu_z * (u_z + w_x) + mu * (u_zz + w_xz);
tauyz_y = mu_y * (v_z + w_y) + mu * (v_yz + w_yy);
tauyz_z = mu_z * (v_z + w_y) + mu * (v_zz + w_yz);

pu_x = p_x * u + p * u_x;
pv_y = p_y * v + p * v_y;
pw_z = p_z * w + p * w_z;
utauxx_x = u_x * tauxx + u * tauxx_x;
vtauxy_x = v_x * tauxy + v * tauxy_x;
wtauxz_x = w_x * tauxz + w * tauxz_x;
utauxy_y = u_y * tauxy + u * tauxy_y;
vtauyy_y = v_y * tauyy + v * tauyy_y;
wtauyz_y = w_y * tauyz + w * tauyz_y;
utauxz_z = u_z * tauxz + u * tauxz_z;
vtauyz_z = v_z * tauyz + v * tauyz_z;
wtauzz_z = w_z * tauzz + w * tauzz_z;

Q_rho  = rho_t  + rhou_x + rhov_y + rhow_z;
Q_rhou = rhou_t + rhouu_x + rhouv_y + rhouw_z + p_x - tauxx_x - tauxy_y - tauxz_z;
Q_rhov = rhov_t + rhouv_x + rhovv_y + rhovw_z + p_y - tauxy_x - tauyy_y - tauyz_z;
Q_rhow = rhow_t + rhouw_x + rhovw_y + rhoww_z + p_z - tauxz_x - tauyz_y - tauzz_z;
Q_rhoe = rhoe_t + rhoue_x + rhove_y + rhowe_z + pu_x + pv_y + pw_z + qx_x + qy_y + qz_z
                - utauxx_x - vtauxy_x - wtauxz_x
                - utauxy_y - vtauyy_y - wtauyz_y
                - utauxz_z - vtauyz_z - wtauzz_z;
