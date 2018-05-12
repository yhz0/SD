\* Source Pyomo model name=Scenario1 *\

min 
_PySP_UserCostObjective:
+4 x1
+2 x2
+10 u1
+10 u2
+0.20000000000000001 v1
+0.20000000000000001 v2
-8 w11
-4 w12
-4 w22
+1 ONE_VAR_CONSTANT

s.t.

r_l_bounds(1)_:
+1 x1
>= 0

r_u_bounds(1)_:
+1 x1
<= 217

r_l_bounds(2)_:
+1 x2
>= 0

r_u_bounds(2)_:
+1 x2
<= 217

c_e_s1_:
-1 x1
+1 v1
+1 w11
+1 w12
= 0

c_e_s2_:
-1 x2
+1 v2
+1 w22
= 0

c_l_bounds(3)_:
+1 v1
>= 0

c_l_bounds(4)_:
+1 v2
>= 0

c_l_bounds(5)_:
+1 u1
>= 0

c_l_bounds(6)_:
+1 u2
>= 0

c_l_bounds(7)_:
+1 w11
>= 0

c_l_bounds(8)_:
+1 w12
>= 0

c_l_bounds(9)_:
+1 w22
>= 0

c_l_d2_lower_:
+1 u2
+1 w12
+1 w22
>= 5.9603195920000003

c_u_d2_upper_:
+1 u2
+1 w12
+1 w22
<= 5.9603195920000003

r_l_d1_:
+1 u1
+1 w11
>= 17.757318649999998

r_u_d1_:
+1 u1
+1 w11
<= 17.757318649999998

c_e_ONE_VAR_CONSTANT: 
ONE_VAR_CONSTANT = 1.0

bounds 
    -inf <= x1 <= +inf
    -inf <= x2 <= +inf
    -inf <= v1 <= +inf
    -inf <= v2 <= +inf
    -inf <= u1 <= +inf
    -inf <= u2 <= +inf
    -inf <= w11 <= +inf
    -inf <= w12 <= +inf
    -inf <= w22 <= +inf
end 
