%% Input Structural Properties
clear 
clc
m = 1;  n = 1;  s = 2;
a = 10; b = 10; c = 5;

%% Nodes Coordinats
n_nodes = (m+1)*(n+1)*(s+1); % Total Number of nodes
Nodes_i = zeros(n_nodes,3); % Coordinates of all nodes
for i = 1:n_nodes
    Nodes_i(i,1) = (rem(i-1,m+1))*a;
    Nodes_i(i,2) = (ceil((rem(i-1,((m+1)*(n+1)))+1)/(m+1))-1)*b;
    Nodes_i(i,3) = (s + 1 - ceil(i/((m+1)*(n+1))))*c;
end

adjacency = zeros(n_nodes,n_nodes);
for i = 1:n_nodes
    for j = i+1:n_nodes
        xA = Nodes_i(i,1);
        yA = Nodes_i(i,2);
        zA = Nodes_i(i,3);
        xB = Nodes_i(j,1);
        yB = Nodes_i(j,2);
        zB = Nodes_i(j,3);
        if (xB-xA == a && yB-yA == 0 && zA-zB == 0 && zA ~= 0) ... 
                || (xB-xA == 0 && yB-yA == b && zA-zB == 0 && zA ~= 0) ...
                || (xB-xA == 0 && yB-yA == 0 && zA-zB == c)
            adjacency(i,j) = 1;
        end
    end
end

%% Model Update after Perturbation
e = (m+1)*(n+1)*s;
delta_Nodes = zeros(n_nodes,3); % Change in elevation of nodes
aa = 1.0;
bb = -1.0;
for i=1:e
    delta_Nodes(i,1) = (bb-aa).*rand(1) + aa;
    delta_Nodes(i,2) = (bb-aa).*rand(1) + aa;
    delta_Nodes(i,3) = (bb-aa).*rand(1) + aa;
end

%% Distributed Algorithm

% Coordinates after perturbation
Nodes_F = Nodes_i + delta_Nodes;
alpha = 0.49; % alpha value

% Difference between final and initial x-, y- and z- coordinates
z = zeros(e,3);
for i=1:e
   z(i,1) = Nodes_F(i,1) - Nodes_i(i,1);
   z(i,2) = Nodes_F(i,2) - Nodes_i(i,2);
   z(i,3) = Nodes_F(i,3) - Nodes_i(i,3);
end

v = (m+1)*(n+1);
nb = zeros(n_nodes,n_nodes);
for k=1:s
    for i=1+(k-1)*v:k*v
        for j=1+(k-1)*v:k*v
            if i==j
                nb(i,j) = 0;
            else
                if adjacency(i,j)==1
                    nb(i,j) = 1;
                    nb(j,i) = 1;
                end
            end
        end
    end
end

nb_v = zeros(n_nodes, n_nodes);
% Interacting column matrix
for i=1:n_nodes
    if (i-v)>0
        nb_v(i,i-v) = 1;
    end
    if (i+v)<n_nodes+1
    nb_v(i,i+v) =1;
    end
end

nb_s = zeros(n_nodes,n_nodes);
% Interacting column matrix
for i=1:n_nodes
    if (i+v)<n_nodes+1
    nb_s(i,i+v) =1;
    end
end

%% Iteration

iteration = 4999;
alpha = zeros(iteration,1);
for k=1:iteration
    alpha(k) = 0.0001*k;
end

%%
tt = zeros(iteration,1);
for k=1:iteration
    t = 0;
    Nodes_f = Nodes_F;
    count = 1;
    while count>0.001
        for i=1:e
            delta_x = 0;
            delta_y = 0;
            delta_z1 = 0;
            delta_z2 = 0;
            for j=1:n_nodes
                delta_x = delta_x + nb_v(i,j)*(Nodes_f(j,1) - Nodes_f(i,1));
                delta_y = delta_y + nb_v(i,j)*(Nodes_f(j,2) - Nodes_f(i,2));
                delta_z1 = delta_z1 + nb(i,j) * (Nodes_f(j,3) - Nodes_f(i,3));
%                 delta_z2 = delta_z2 + nb_s(i,j)*(c - abs(Nodes_f(i,3) - Nodes_f(j,3)));
            end
            z(i,1) = delta_x;
            z(i,2) = delta_y;
            z(i,3) = delta_z1 + delta_z2;
        end
        count = max(max(abs(z)));
        t = t + 1;
        Nodes_f(1:e,:) = Nodes_f(1:e,:) + alpha(k)*z;
    end
    tt(k) = t;
end

C = [tt, alpha]

