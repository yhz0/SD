\* Source Pyomo model name=Scenario1 *\

min 
x1:
+4 x16
+2 x17
+0.20000000000000001 x18
+0.20000000000000001 x19
+10 x20
+10 x21
-8 x22
-4 x23
-4 x24
+1 ONE_VAR_CONSTANT

s.t.

r_l_x7_:
+1 x16
>= 0

r_u_x7_:
+1 x16
<= 217

r_l_x8_:
+1 x17
>= 0

r_u_x8_:
+1 x17
<= 217

c_e_x2_:
-1 x17
+1 x19
+1 x24
= 0

c_e_x5_:
-1 x16
+1 x18
+1 x22
+1 x23
= 0

c_l_x10_:
+1 x19
>= 0

c_l_x11_:
+1 x20
>= 0

c_l_x12_:
+1 x21
>= 0

c_l_x13_:
+1 x22
>= 0

c_l_x14_:
+1 x23
>= 0

c_l_x15_:
+1 x24
>= 0

c_l_x3_:
+1 x21
+1 x23
+1 x24
>= 5.9603195920000003

c_l_x9_:
+1 x18
>= 0

c_u_x4_:
+1 x21
+1 x23
+1 x24
<= 5.9603195920000003

r_l_x6_:
+1 x20
+1 x22
>= 17.757318649999998

r_u_x6_:
+1 x20
+1 x22
<= 17.757318649999998

c_e_ONE_VAR_CONSTANT: 
ONE_VAR_CONSTANT = 1.0

bounds 
    -inf <= x16 <= +inf
    -inf <= x17 <= +inf
    -inf <= x18 <= +inf
    -inf <= x19 <= +inf
    -inf <= x20 <= +inf
    -inf <= x21 <= +inf
    -inf <= x22 <= +inf
    -inf <= x23 <= +inf
    -inf <= x24 <= +inf
end 
