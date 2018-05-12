\* Source Pyomo model name=Scenario1 *\

min 
x1:
+4 x6
+2 x7
+10 x10
+10 x11
-8 x12
-4 x13
-4 x14
+0.20000000000000001 x8
+0.20000000000000001 x9

s.t.

c_e_x2_:
+1 x10
+1 x12
= 17.757318649999998

c_e_x3_:
+1 x11
+1 x13
+1 x14
= 5.9603195920000003

c_e_x4_:
-1 x6
+1 x12
+1 x13
+1 x8
= 0

c_e_x5_:
-1 x7
+1 x14
+1 x9
= 0

c_e_ONE_VAR_CONSTANT: 
ONE_VAR_CONSTANT = 1.0

bounds 
   0 <= x6 <= 217
   0 <= x7 <= 217
   0 <= x8 <= +inf
   0 <= x9 <= +inf
   0 <= x10 <= +inf
   0 <= x11 <= +inf
   0 <= x12 <= +inf
   0 <= x13 <= +inf
   0 <= x14 <= +inf
end 
