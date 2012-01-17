# 
# 
# 
# 
# Blotnner Viscosity Models for N and N2 --------------------------------------------------------------------------------------------

mu_N := (1/10)*exp((A_N*ln(T)+B_N)*ln(T)+C_N); mu_N2 := (1/10)*exp((A_N2*ln(T)+B_N2)*ln(T)+C_N2);


# Wilke mixing function -------------------------------------------------------------------------------------------------------------
;

c := [c_N, c_N2]; M := [M_N, M_N2]; MU := [mu_N, mu_N2];
Wilke_mixing := proc (VectorProperty, c, M) local chi, Mtot, s, r, temp, num, den, phi, Property_mix; Mtot := 0; for s to 2 do chi[s] := c[s]/M[s]; Mtot := Mtot+chi[s] end do; Mtot := 1/Mtot; for s to 2 do chi[s] := chi[s]*Mtot end do; for s to 2 do phi[s] := 0; for r to 2 do temp := 1+sqrt(MU[s]/MU[r])*sqrt(sqrt(M[r]/M[s])); num := chi[r]*temp*temp; den := sqrt(8+8*M[s]/M[r]); phi[s] := phi[s]+num/den end do end do; Property_mix := 0.; for s to 2 do Property_mix := Property_mix+VectorProperty[s]*chi[s]/phi[s] end do end proc;


# Applying Wilke_mixing to the Blotneer viscosity models of N and N2: 
mu_mix := Wilke_mixing(MU, c, M);

# 
