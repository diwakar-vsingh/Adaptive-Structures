%% Input Structural Properties
clear 
clc
m = 1;  n = 1;  s = 2;
a = 10; b = 10; c = 5;

%% Nodes Coordinats
n_nodes = (m+1)*(n+1)*(s+1); % Total Number of nodes
iteration = 10000;
Nodes_i = zeros(n_nodes,3,iteration); % Coordinates of all nodes
for k=1:iteration
    for i = 1:n_nodes
        Nodes_i(i,1,k) = (rem(i-1,m+1))*a;
        Nodes_i(i,2,k) = (ceil((rem(i-1,((m+1)*(n+1)))+1)/(m+1))-1)*b;
        Nodes_i(i,3,k) = (s + 1 - ceil(i/((m+1)*(n+1))))*c;
    end
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

% This returns random value in range [-1,0.5] 
aa = 1.0;
bb = -1.0;
for k=1:iteration
   for i=1:e
        delta_Nodes(i,1,k) = (bb-aa).*rand(1) + aa;
        delta_Nodes(i,2,k) = (bb-aa).*rand(1) + aa;
        delta_Nodes(i,3,k) = (bb-aa).*rand(1) + aa;
    end 
end

%% Distributed Algorithm

% Coordinates after perturbation
Nodes_F = Nodes_i + delta_Nodes;
alpha = 0.35; % alpha value

% Difference between final and initial x-, y- and z- coordinates
z = zeros(e,3,iteration);
for k=1:iteration
    for i=1:e
        z(i,1,k) = Nodes_F(i,1,k) - Nodes_i(i,1,k);
        z(i,2,k) = Nodes_F(i,2,k) - Nodes_i(i,2,k);
        z(i,3,k) = Nodes_F(i,3,k) - Nodes_i(i,3,k);
    end
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

Nodes_f = Nodes_F;
X_t(1,:,:) = Nodes_i(:,1,:);
Y_t(1,:,:) = Nodes_i(:,2,:);
Z_t(1,:,:) = Nodes_i(:,3,:);
X_t(2,:,:) = Nodes_f(:,1,:);
Y_t(2,:,:) = Nodes_f(:,2,:);
Z_t(2,:,:) = Nodes_f(:,3,:);

%%
tt = zeros(iteration,1);
for k=1:iteration
    check = 1;
    t = 0;
    while check>0.001
        for i=1:e
            delta_x = 0;
            delta_y = 0;
            delta_z1 = 0;
            delta_z2 = 0;
            for j=1:n_nodes
                delta_x = delta_x + nb_v(i,j) * (Nodes_f(j,1,k) - Nodes_f(i,1,k));
                delta_y = delta_y + nb_v(i,j) * (Nodes_f(j,2,k) - Nodes_f(i,2,k));
                delta_z1 = delta_z1 + nb(i,j) * (Nodes_f(j,3,k) - Nodes_f(i,3,k));
%                 delta_z2 = delta_z2 + nb_s(i,j)*(c - abs(Nodes_f(i,3,k) - Nodes_f(j,3,k)));
            end
                z(i,1,k) = delta_x;
                z(i,2,k) = delta_y;
                z(i,3,k) = delta_z1 + delta_z2;
        end
        c1 = z(:,:,k);
        check = max(max(abs(c1)));
        t = t + 1;
        Nodes_f(1:e,:,k) = Nodes_f(1:e,:,k) + alpha*z(:,:,k);
%         X_t(t(k)+2,:,k) = Nodes_f(:,1,k);
%         Y_t(t(k)+2,:,k) = Nodes_f(:,2,k);
%         Z_t(t(k)+2,:,k) = Nodes_f(:,3,k);
    end
    tt(k) = t;
end
[min(t) max(t)]
