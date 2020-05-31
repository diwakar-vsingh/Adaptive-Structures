%% Input Structural Properties
clear 
clc
m = 1; % Number of bays in x-direction
n = 1; % Number of bays in y-direction
s = 2; % Number of storeys 
a = 10; % Dimension of bay along x-direction
b = 10; % Dimension of bay along y-direction
c = 5; % Height of storey

%% Nodes Coordinates
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

%% Figure

figure
% Drawing line plot in 3D
for i=1:n_nodes
    for j=1:n_nodes
        line=adjacency(i,j);
        if line ~= 0 
            plot3([Nodes_i(i,1),Nodes_i(j,1)],[Nodes_i(i,2),Nodes_i(j,2)],[Nodes_i(i,3),Nodes_i(j,3)],'k');
            hold all;
        end
    end    
end
% % Numbering all nodes
% for i=1:n_nodes
%      text(Nodes_i(i,1),Nodes_i(i,2),Nodes_i(i,3),['\bf',num2str(i)],'Color',[0 0.0 1.0],'FontSize',14);
%      hold all;
% end

xlabel('x-axis')
ylabel('y-axis')
title('3D Frame Structure')
clc

%% Model Update after Perturbation
e = (m+1)*(n+1)*s;
delta_Nodes = zeros(n_nodes,3); % Change in elevation of nodes

% dlg_title = 'Perturbation';
% prompt = cell(1,e);
% defaultans = cell(1,e);
% for i=1:e
%     prompt(i) = {string(['Lateral Displacement in x-direction at Node ' num2str(i) ':'])};
%     defaultans(i) = {'0'};
% end
% x = inputdlg(prompt, dlg_title, 1, defaultans); 
% 
% for i=1:e
%     prompt(i) = {string(['Lateral Displacement in y-direction at Node ' num2str(i) ':'])};
%     defaultans(i) = {'0'};
% end
% y = inputdlg(prompt, dlg_title, 1, defaultans);
% 
% for i=1:e
%     prompt(i) = {string(['Vertical Displacement at Node ' num2str(i) ':'])};
%     defaultans(i) = {'0'};
% end
% z = inputdlg(prompt, dlg_title, 1, defaultans); 
% 
% % Model Update
% for i=1:e
%     delta_Nodes(i,1) =  str2double(x{i,:}); 
%     delta_Nodes(i,2) =  str2double(y{i,:}); 
%     delta_Nodes(i,3) =  str2double(z{i,:}); 
% end

% This returns random value in range [-1,0.5] for lateral displacement and
% between [-0.5, 0.25]
aa = 1.0;
bb = -1.0;
for i=1:e
    delta_Nodes(i,1) = (bb-aa).*rand(1) + aa;
    delta_Nodes(i,2) = (bb-aa).*rand(1) + aa;
    delta_Nodes(i,3) = (bb-aa).*rand(1) + aa;
end

%% Distributed Algorithm

% Nodes Coordinates after perturbation
Nodes_F = Nodes_i + delta_Nodes;
alpha = 0.35; % alpha value
% Difference between final and initial x-, y- and z- coordinates
z = zeros(e,3);
for i=1:e
   z(i,1) = Nodes_F(i,1) - Nodes_i(i,1);
   z(i,2) = Nodes_F(i,2) - Nodes_i(i,2);
   z(i,3) = Nodes_F(i,3) - Nodes_i(i,3);
end

v = (m+1)*(n+1);
nb = zeros(n_nodes,n_nodes); % Neighboring nodes for controlling slab tilt
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

nb_v = zeros(n_nodes, n_nodes); % Neighboring nodes for twist and/or drift
% Interacting column matrix
for i=1:n_nodes
    if (i-v)>0
        nb_v(i,i-v) = 1;
    end
    if (i+v)<n_nodes+1
    nb_v(i,i+v) =1;
    end
end

nb_s = zeros(n_nodes,n_nodes); % Neighboring nodes for maintaining initial storey height
% Interacting column matrix
for i=1:n_nodes
    if (i+v)<n_nodes+1
    nb_s(i,i+v) =1;
    end
end

% Drawing line plot in 3D
for i=1:n_nodes
    for j=1:n_nodes
        line=adjacency(i,j);
        if line ~= 0 
            plot3([Nodes_F(i,1),Nodes_F(j,1)],[Nodes_F(i,2),Nodes_F(j,2)],[Nodes_F(i,3),Nodes_F(j,3)],'--r');
            hold all;
        end
    end    
end
hold off
axis off

%% Iteration
t = 0;
Nodes_f = Nodes_i + delta_Nodes;
X_t(1,:) = Nodes_i(:,1);
Y_t(1,:) = Nodes_i(:,2);
Z_t(1,:) = Nodes_i(:,3);
X_t(2,:) = Nodes_f(:,1);
Y_t(2,:) = Nodes_f(:,2);
Z_t(2,:) = Nodes_f(:,3);
while max(max(abs(z)))>0.001
    for i=1:e
        delta_x = 0;
        delta_y = 0;
        delta_z1 = 0;
        delta_z2 = 0;
        for j=1:n_nodes
            delta_x = delta_x + nb_v(i,j)*(Nodes_f(j,1) - Nodes_f(i,1));
            delta_y = delta_y + nb_v(i,j)*(Nodes_f(j,2) - Nodes_f(i,2));
            delta_z1 = delta_z1 + nb(i,j) * (Nodes_f(j,3) - Nodes_f(i,3));
%             delta_z2 = delta_z2 + nb_s(i,j)*(c - abs(Nodes_f(i,3) - Nodes_f(j,3)));
        end
        z(i,1) = delta_x;
        z(i,2) = delta_y;
        z(i,3) = delta_z1 + delta_z2;
    end
    t = t+1;
    Nodes_f(1:e,:) = Nodes_f(1:e,:) + alpha*z;
    X_t(t+2,:) = Nodes_f(:,1);
    Y_t(t+2,:) = Nodes_f(:,2);
    Z_t(t+2,:) = Nodes_f(:,3);
end
C = zeros(t+2,3,n_nodes);
for i=1:n_nodes
    C(:,:,i) = [X_t(:,i), Y_t(:,i), Z_t(:,i)];
end
% Drawing line plot in 3D
figure
for i=1:n_nodes
    for j=1:n_nodes
        line=adjacency(i,j);
        if line ~= 0 
            plot3([Nodes_i(i,1),Nodes_i(j,1)],[Nodes_i(i,2),Nodes_i(j,2)],[Nodes_i(i,3),Nodes_i(j,3)],'k');
            hold all;
            plot3([Nodes_F(i,1),Nodes_F(j,1)],[Nodes_F(i,2),Nodes_F(j,2)],[Nodes_F(i,3),Nodes_F(j,3)],'--r');
            hold all;
            plot3([Nodes_f(i,1),Nodes_f(j,1)],[Nodes_f(i,2),Nodes_f(j,2)],[Nodes_f(i,3),Nodes_f(j,3)],'-.k');
            hold all;
        end
    end    
end
xlabel('x-axis')
ylabel('y-axis')
title('Final Configuration')
clc

%%
Comp = [Nodes_i, Nodes_f];
['Number of iterations: ' num2str(t)]
x = t-1;
figure
subplot(3,2,1)
plot(X_t(:,[1,3,5,7]),'DisplayName','X_t(:,[1,3,5,7])')
hold on
plot([0,x],[0,0], ':')
xlabel('Iteration')
ylabel('x(t)')

subplot(3,2,2)
plot(X_t(:,[2,4,6,8]),'DisplayName','X_t(:,[2,4,6,8])')
hold on
plot([0,x],[10,10], ':')
xlabel('Iteration')
ylabel('x(t)')

subplot(3,2,3)
plot(Y_t(:,[1:2,5:6]),'DisplayName','Y_t(:,[1:2,5:6])')
hold on
plot([0,x],[0,0], ':')
xlabel('Iteration')
ylabel('y(t)')

subplot(3,2,4)
plot(Y_t(:,[3:4,7:8]),'DisplayName','Y_t(:,[3:4,7:8])')
hold on
plot([0,x],[10,10], ':')
xlabel('Iteration')
ylabel('y(t)')

subplot(3,2,5)
plot(Z_t(:,1:4),'DisplayName','Z_t(:,1:4)')
hold on
plot([0,x],[Z_t(x,4),Z_t(45,4)], ':')
xlabel('Iteration')
ylabel('z(t)')

subplot(3,2,6)
plot(Z_t(:,5:8),'DisplayName','Z_t(:,5:8)')
hold on
plot([0,x],[Z_t(x,8),Z_t(45,8)], ':')
xlabel('Iteration')
ylabel('z(t)')
title('Convergence of x-,y-, and z-coordinate')

