clc
clear
H = 10;
B = 15;
L = 20;
C_initial = [0 H 0; L H 0; 0 H B; L H B];
Y = [0 -1 0; 0 -0.5 0; 0 0.75 0; 0 -0.5 0; ];
delta_y = 0.2;
C_final = C_initial + delta_y*Y;
y_i = C_final(1,2); 
y_j = C_final(2,2);
y_k = C_final(3,2); 
y_l = C_final(4,2);
alpha = 0.3;
z_i = C_final(1,2) - C_initial(1,2);
z_j = C_final(2,2) - C_initial(2,2);
z_k = C_final(3,2) - C_initial(3,2);
z_l = C_final(4,2) - C_initial(4,2);
t = 0;

while abs(z_i)>0.001 || abs(z_j)>0.001 || abs(z_k)>0.001 || abs(z_l)>0.001
    z_i = (y_j-y_i) + (y_k-y_i) + (y_l-y_i);
    z_j = (y_i-y_j) + (y_k-y_j) + (y_l-y_j);
    z_k = (y_i-y_k) + (y_j-y_k) + (y_l-y_k);
    z_l = (y_i-y_l) + (y_j-y_l) + (y_k-y_l);
    y_i = y_i + alpha*z_i
    y_j = y_j + alpha*z_j
    y_k = y_k + alpha*z_k
    y_l = y_l + alpha*z_l
    t = t+1
end
