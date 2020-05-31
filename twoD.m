clc
clear
H = 10;
L = 20;
delta_y = 0.5;
C_initial = [0 H; L H];
Y = [0 1; 0 -2];
C_final = C_initial + delta_y*Y;
y_i = C_final(1,2); 
y_j = C_final(2,2);
alpha = 0.5;
z_i = C_final(1,2) - C_initial(1,2);
z_j = C_final(2,2) - C_initial(2,2);
t = 0;

while abs(z_i)>0.0001 && abs(z_j)>0.0001
    z_i = alpha*(y_j-y_i);
    z_j = alpha*(y_i-y_j);
    y_i = y_i + z_i;
    y_j = y_j + z_j;
    t = t+1;
end
