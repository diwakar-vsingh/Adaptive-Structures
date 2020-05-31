%% Input Structural Properties
clear 
clc
% y=0;
% while y == 0
%     prompt = {'Number of bays along x-direction: ','Number of bays along y-direction: ',...
%         'Number of storey: ', 'Dimension of one bay along length: ', 'Dimension of one bay along breadth: ',...
%         'Dimension of one storey: '};
%     dlg_title = 'Structural Properties';
%     defaultans = {'1','1','1','10','5','10'};
%     x = inputdlg(prompt,dlg_title,[1 50; 1 50; 1 50; 1 50; 1 50; 1 50], defaultans);
%     m = str2double(x{1,:}); 
%     n = str2double(x{2,:}); 
%     s = str2double(x{3,:});
%     a = str2double(x{4,:}); 
%     b = str2double(x{5,:}); 
%     c = str2double(x{6,:});
%     if s>0 && m>=0 && n>=0 && a>0 && b>0 && c>0 && rem(s,1)==0 && rem(m,1)==0 && rem(n,1)==0
%         y = 1;
%     elseif s<=0 || a<=0 || b<=0 || c<=0
%         uiwait(warndlg('Number of Storey or member lengths cannot be zero or negative','!! Error !!','modal'));
%         y = 0;
%     elseif rem(m,1)~=0 || rem(n,1)~=0
%         uiwait(warndlg('Number of bays along length or width cannot by fractional','!! Error !!','modal'));
%         y = 0;
%     elseif rem(s,1)~=0 
%         uiwait(warndlg('Number of Storey cannot by fractional','!! Error !!','modal'));
%         y = 0;
%     elseif m<0 || n<0
%         uiwait(warndlg('Number of bays along length or width cannot be negative','!! Error !!','modal'));
%         y = 0;
%     else
%         uiwait(warndlg('Please input numbers','!! Error !!','modal'));
%         y = 0;        
%     end
% end
m = 1; n = 1; s = 2;
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

%% Figure

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
% Numbering all nodes
for i=1:n_nodes
     text(Nodes_i(i,1),Nodes_i(i,2),Nodes_i(i,3),['\bf',num2str(i)],'Color',[0 0.0 1.0],'FontSize',14);
     hold all;
end
% Drawing fixed supports
for i=0:m
    for j=0:n
        text(i*a-0.05*a,j*b,0.07*c,'__','Color',[0 0.0 0.0],'FontSize',10);
        text(i*a-0.05*a,j*b,0.02*c,',,,,','Color',[0 0.0 0.0],'FontSize',10);
        hold all;
    end
end
xlabel('x-axis')
ylabel('y-axis')
title('3D Frame Structure')


%% Degree of freedom associated with each nodes
dof=zeros(n_nodes,6);
for i=1:(m+1)*(n+1)*(s+1)
    for j=1:6
        dof(i,j)=6*i-6+j;
    end
end
        
%% Association matrix
member_nos = 0;
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
            member_nos = member_nos + 1;
            if xB-xA == a && yB-yA == 0 && zA-zB == 0
                L(member_nos) = a;
                t(member_nos) = 1; % Element in x-direction
            elseif xB-xA == 0 && yB-yA == b && zA-zB == 0
                L(member_nos) = b;
                t(member_nos)= 2 ; % Element in y-direction
            elseif xB-xA == 0 && yB-yA == 0 && zA-zB == c 
                L(member_nos) = c;
                t(member_nos) = 3; % Element in z-direction
            end
            elements(member_nos,1) = i;
            elements(member_nos,2) = j;
            asscn(member_nos,:) = [dof(i,1) dof(i,2) dof(i,3) dof(i,4) dof(i,5) dof(i,6)...
                dof(j,1) dof(j,2) dof(j,3) dof(j,4) dof(j,5) dof(j,6)];            
        end
    end
end

%% Model Update after Perturbation
y = 0;
e = (m+1)*(n+1)*s;
delta_Nodes = zeros(n_nodes,3); % Change in elevation of nodes
dlg_title = 'Perturbation';
prompt = cell(1,e);
defaultans = cell(1,e);

for i=1:e
    prompt(i) = {string(['Change in Elevation at Node ' num2str(i) ':'])};
    defaultans(i) = {'0'};
end
x = inputdlg(prompt,dlg_title,1,defaultans); % Input of perturbation

% Model Update
for i=1:e
    delta_Nodes(i,3) =  str2double(x{i,:}); 
end

for i=1:e
    a = 0.5;
    b = -1.0;
    delta_Nodes(i,3) = (b-a).*rand(1) + a;
end
% for i=1:e
%     if rem(i,2)==0
%         delta_Nodes(i,3) = -1;
%     else
%         delta_Nodes(i,3) = 1;
%     end
% end

%% Distributed Algorithm
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

nb_s = zeros(n_nodes,n_nodes);
% Interacting column matrix
for i=1:n_nodes
    if (i-v)>0
        nb_s(i,i-v) = 1;
    end
    if (i+v)<n_nodes+1
    nb_s(i,i+v) =1;
    end
end

%%
clc
Nodes_f = Nodes_i + delta_Nodes; % Coordinates after perturbation
alpha = 0.2; % alpha value
% Drawing line plot in 3D
for i=1:n_nodes
    for j=1:n_nodes
        line=adjacency(i,j);
        if line ~= 0 
            plot3([Nodes_f(i,1),Nodes_f(j,1)],[Nodes_f(i,2),Nodes_f(j,2)],[Nodes_f(i,3),Nodes_f(j,3)],'--r');
            hold all;
        end
    end    
end

%% Iteration
z = zeros(e,1);
for i=1:e
   z(i) = Nodes_f(i,3) - Nodes_i(i,3); 
end
t = 0;
while max(abs(z))>0.001
    for i=1:e
        delta1 = 0;
        delta2 = 0;
        for j=1:n_nodes
            delta1 = delta1 + nb(i,j)*(Nodes_f(j,3) - Nodes_f(i,3));
        end
        z(i) = delta1 + delta2;
    end
%     Comp(t+1,1) = Nodes_i(1,3);
%     Comp(t+1,2) = Nodes_f(1,3);
    Nodes_f(1:e,3) = Nodes_f(1:e,3) + alpha*z;
    t = t+1;
end

%% Drawing line plot in 3D
for i=1:n_nodes
    for j=1:n_nodes
        line=adjacency(i,j);
        if line ~= 0 
            plot3([Nodes_f(i,1),Nodes_f(j,1)],[Nodes_f(i,2),Nodes_f(j,2)],[Nodes_f(i,3),Nodes_f(j,3)],'-.k');
            hold all;
        end
    end    
end
clc
Cc = [Nodes_i(1:e,3), Nodes_f(1:e,3)]
t




