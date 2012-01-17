# 

# Relation of molar masses of N and N2:
NULL;
NULL;
NULL;


# Mass fraction of Nitrogen atom and molecule - c_N and c_N2
NULL;

# Translational/rotational specific heat :--------------------------
Cv_tr_N := (3/2)*R_N; Cv_tr_N2 := (5/2)*R_N2;

# Translational specific heat :-------------------------------------
Cv_trans_N := (3/2)*R_N; Cv_trans_N2 := (3/2)*R_N2;

# Rotational specific heat :----------------------------------------
Cv_rot_N := Cv_tr_N-Cv_trans_N; Cv_rot_N2 := Cv_tr_N2-Cv_trans_N2;

# Vibrational specific heat : --------------------------------------
e_vib_N := 0; Cv_vib_N := 0; e_vib_N2 := R_N2*theta_v_N2/(exp(theta_v_N2/Tv)-1); Cv_vib_N2 := R_N2*theta_v_N2^2*exp(theta_v_N2/Tv)/((exp(theta_v_N2/Tv)-1)^2*Tv^2);
# Electronic specific heat :---------------------------------------------------
e_elec_N := R_N*(sum(theta_e_N[i]*g_N[i]*exp(-theta_e_N[i]/Tv), i = 0 .. energy_level_N))/(sum(g_N[i]*exp(-theta_e_N[i]/Tv), i = 0 .. energy_level_N)); Cv_elec_N := R_N*(sum(theta_e_N[i]^2*g_N[i]*exp(-theta_e_N[i]/Tv)/Tv^2, i = 0 .. energy_level_N))/(sum(g_N[i]*exp(-theta_e_N[i]/Tv), i = 0 .. energy_level_N))-R_N*(sum(theta_e_N[i]*g_N[i]*exp(-theta_e_N[i]/Tv), i = 0 .. energy_level_N))*(sum(theta_e_N[i]*g_N[i]*exp(-theta_e_N[i]/Tv)/Tv^2, i = 0 .. energy_level_N))/(sum(g_N[i]*exp(-theta_e_N[i]/Tv), i = 0 .. energy_level_N))^2; e_elec_N2 := R_N2*(sum(theta_e_N2[i]*g_N2[i]*exp(-theta_e_N2[i]/Tv), i = 0 .. energy_level_N2))/(sum(g_N2[i]*exp(-theta_e_N2[i]/Tv), i = 0 .. energy_level_N2)); Cv_elec_N2 := R_N2*(sum(theta_e_N2[i]^2*g_N2[i]*exp(-theta_e_N2[i]/Tv)/Tv^2, i = 0 .. energy_level_N2))/(sum(g_N2[i]*exp(-theta_e_N2[i]/Tv), i = 0 .. energy_level_N2))-R_N2*(sum(theta_e_N2[i]*g_N2[i]*exp(-theta_e_N2[i]/Tv), i = 0 .. energy_level_N2))*(sum(theta_e_N2[i]*g_N2[i]*exp(-theta_e_N2[i]/Tv)/Tv^2, i = 0 .. energy_level_N2))/(sum(g_N2[i]*exp(-theta_e_N2[i]/Tv), i = 0 .. energy_level_N2))^2;
# 
# Four modes of conductivity for N --------------------------------------------------------------------------------------------------
# 
k_trans_N := (5/2)*mu_N*Cv_trans_N;
k_rot_N := mu_N*Cv_rot_N;
k_vib_N := mu_N*Cv_vib_N;
k_elec_N := mu_N*Cv_elec_N;

# Four modes of conductivity for N --------------------------------------------------------------------------------------------------
k_trans_N2 := (5/2)*mu_N2*Cv_trans_N2;
k_rot_N2 := mu_N2*Cv_rot_N2;
k_vib_N2 := mu_N2*Cv_vib_N2;
k_elec_N2 := mu_N2*Cv_elec_N2;

# Wilke mixing function -------------------------------------------------------------------------------------------------------------
;

c := [c_N, c_N2]; M := [M_N, M_N2]; Kappa_e := [k_elec_N, k_elec_N2]; Kappa_v := [k_vib_N, k_vib_N2]; Kappa_t := [k_trans_N, k_trans_N2]; Kappa_r := [k_rot_N, k_rot_N2]; MU := [mu_N, mu_N2];
# It reproduces the procedure "void Transport::MixtureTransport::init_wilke_mixing () const" arround line 251 of file: ~/FINS/trunk/physics/properties/transport.C

Wilke_mixing := proc (VectorProperty, c, M) local chi, Mtot, s, r, temp, num, den, phi, Property_mix; Mtot := 0; for s to 2 do chi[s] := c[s]/M[s]; Mtot := Mtot+chi[s] end do; Mtot := 1/Mtot; for s to 2 do chi[s] := chi[s]*Mtot end do; for s to 2 do phi[s] := 0; for r to 2 do temp := 1+sqrt(MU[s]/MU[r])*sqrt(sqrt(M[r]/M[s])); num := chi[r]*temp*temp; den := sqrt(8+8*M[s]/M[r]); phi[s] := phi[s]+num/den end do end do; Property_mix := 0.; for s to 2 do Property_mix := Property_mix+VectorProperty[s]*chi[s]/phi[s] end do end proc;
# Note that Wilke mixing is linear in k, so Wilke_mixing(Kappa_e+Kappa_v, c, M) == Wilke_mixing(Kappa_e, c, M)+Wilke_mixing(Kappa_v, c, M)
# I decided to calculate the Wilke mixing for N and N2 for of each mode. And then adding it where necessary/convenient.


kv_mix := Wilke_mixing(Kappa_v, c, M);
ke_mix := Wilke_mixing(Kappa_e, c, M);
kt_mix := Wilke_mixing(Kappa_t, c, M);
kr_mix := Wilke_mixing(Kappa_r, c, M);

k_tr_mix := kt_mix+kr_mix;
k_ev_mix := ke_mix+kv_mix;

k_mix := k_tr_mix+k_ev_mix;


NULL;
Kappa := Kappa_t+Kappa_r+Kappa_v+Kappa_e;
k_mix_sum := Wilke_mixing(Kappa, c, M);
Res := simplify(k_mix-k_mix_sum, symbol);
Res := simplify(Res, size);

