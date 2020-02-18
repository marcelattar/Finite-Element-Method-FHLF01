%% FEM projekt, loads all the numerical data aswell as the mesh.
clear all;

load t.mat
load e.mat
load p.mat

p=p/1000; % Converts to meter.

nelm=length(t(1,:));
edof(:,1)=1:nelm;
edof(:,2:4)=t(1:3,:)'; 
coord=p' ;
ndof=max(max(t(1:3,:))) ;
[Ex,Ey]=coordxtr(edof,coord,(1:ndof)',3);

%eldraw2(Ex,Ey,[1,4,1]) 
%grid on

% Numerical data -----------------------
Esmd = 105e9;           Epcb = 105e9;           Esol = 50e9;                    T0 = 30;
vsmd = 0.118;          vpcb = 0.136;           vsol = 0.36;                    ac = 40;
ksmd = 0.29;            kpcb = 1.059;           ksol = 66.8;                    qel = 9e3;
rhosmd = 1850;          rhopcb = 1850;          rhosol = 7265;                  Tinf = 20;
csmd = 950;             cpcb = 950;             csol = 210;
alphasmd = 1.2e-5;      alphapcb = 2e-5;        alphasol = 1.2e-5;
% -------------------------
%% PART 1: Calcualtes the K and C matrix

[rowsEx,~] = size(Ex);
K = zeros(ndof);
C = zeros(ndof);

for i=1:rowsEx
    ex = Ex(i,:);
    ey = Ey(i,:);
    ep = 1;                                 % element thickness 
    if t(4,i)==1                            % Solder (area 1)
        D = ksol*eye(2);
        x = rhosol*csol;                    % x: material paramater for C-matrix
    end
    if t(4,i)==2                            % PCB (area 2) 
        D = kpcb*eye(2);
        x = rhopcb*cpcb;
    end
    if t(4,i)== 3                           % SMD (area 3)
        D = ksmd*eye(2);  
        x = rhosmd*csmd;
    end
    edofrow = edof(i,:);
    Ke = flw2te(ex,ey,ep,D);
    Ce=plantml(ex,ey,x);
    
    K = assem(edofrow,K,Ke);
    C = assem(edofrow,C,Ce);
end
%% Test for isMember
[~,colse] = size(e);
 counterL1 = 0;
 counterL2 = 0;
 counterL3 = 0;
 counterL4 = 0;
 counterL5 = 0;
 counterL6 = 0;
 for i=1:colse
     if isMemberL1(i,e,p)
         counterL1= counterL1 + 1;
     end
     if isMemberL2(i,e,p)
        counterL2 = counterL2 + 1; 
     end
     if isMemberL3(i,e,p)
        counterL3 = counterL3 + 1; 
     end
     if isMemberL4(i,e,p)
        counterL4 = counterL4 + 1; 
     end
     if isMemberL5(i,e,p)
        counterL5 = counterL5 + 1; 
     end
     if isMemberL6(i,e,p)
        counterL6 = counterL6 + 1; 
     end
end
counterL1
counterL2
counterL3
counterL4
counterL5
counterL6
%% Calcualtes the new K-matrx and f-matrix
[~,colse] = size(e);
f = zeros([ndof 1]);
Kc = zeros(ndof);


for i=1:colse
   l = abs(e(3,colse)-e(4,colse))*10^(-3);                          % length of element (converts to meter)
   if isMemberL1(i,e,p)                                             % Test of which boundary the edge is in
       f(e(1,i)) = f(e(1,i)) + qel*l/2;
       f(e(2,i)) = f(e(2,i)) + qel*l/2;
   
   elseif isMemberL2(i,e,p) || isMemberL3(i,e,p);                   % Test of which boundary the edge is in
       f(e(1,i)) = f(e(1,i)) + ac*Tinf*l/2;
       f(e(2,i)) = f(e(2,i)) + ac*Tinf*l/2;
       
       Kc(e(1,i),e(1,i)) = Kc(e(1,i),e(1,i)) + ac*2*l/6;            % The diagonal elements.
       Kc(e(2,i),e(2,i)) = Kc(e(2,i),e(2,i)) + ac*2*l/6;
       
       Kc(e(1,i),e(2,i)) = Kc(e(1,i),e(2,i)) + ac*l/6;              % The non-diagonal elements.
       Kc(e(2,i),e(1,i)) = Kc(e(2,i),e(1,i)) + ac*l/6;

   end
end

Knew = K + Kc;
a=solveq(Knew,f);

%% Plots the stationary solution
ed = extract(edof,a);
%ed = [ed;ed];                                               % This is to mirror the solution
%Ex = [Ex;-Ex];                                              % This is to mirror the solution
%Ey = [Ey;Ey];                                               % This is to mirror the solution

fill(Ex',Ey',ed','EdgeColor','none')                        % Plots the temperature
xlabel('width [meter]', 'fontweight','bold','fontsize',15);
ylabel('height [meter]', 'fontweight','bold','fontsize',15);
title('The stationary temperature [Celsius]','fontsize',16);
set(gca,'fontsize',12);
colorbar;

%% Plots the transient solution
aold = T0*ones([ndof 1]); 
tf = 70;                                                     % For how long you want to go
N = 300;                                                     % Number of steps

a = eulerstep(aold, f, C, Knew, N, tf);

ed = extract(edof,a);
%ed = [ed;ed];                                               % This is to mirror the solution
%Ex = [Ex;-Ex];                                              % This is to mirror the solution
%Ey = [Ey;Ey];                                               % This is to mirror the solution

fill(Ex',Ey',ed','EdgeColor','none')                        % Plots the temperature
xlabel('width [meter]', 'fontweight','bold','fontsize',15);
ylabel('height [meter]', 'fontweight','bold','fontsize',15);
title('The temperature after 70 seconds [Celsius]','fontsize',16);
set(gca,'fontsize',12);
colorbar;

%% PART 2: Compute new Edof-matrix

nedof = [edof(:,1), zeros([nelm 1]), 2*edof(:,2), zeros([nelm 1]), 2*edof(:,3), zeros([nelm 1]), 2*edof(:,4)];  

% This is our new Edof matrix

for i = 1:nelm
   for s = 1:3 
        nedof(:,2*s) = nedof(:,2*s+1)-1;
   end
end

ndof = ndof*2;                                          % The new degree of freedom is 2x bigger

%% Computes the K and f matrix
ptype = 2;                                              % Because of plain strain problem
K = zeros(ndof);
f = zeros([ndof 1]);

for i =1:nelm
    ep = [ptype 1];
    ex = Ex(i,:);
    ey = Ey(i,:);
    
    DeltaTp1 = a(t(1,i))-T0;                            % The temperature difference in each node
    DeltaTp2 = a(t(2,i))-T0;
    DeltaTp3 = a(t(3,i))-T0;
    DeltaT = (DeltaTp1 + DeltaTp2 + DeltaTp3)/3;        % This is the mean temperatuer difference for an element
    
    if t(4,i)==1                                        % Solder
        D = hooke(ptype,Esol,vsol);
        Deps0 = Esol*alphasol*DeltaT*[1;1:1;0]./(1-2*vsol);
    end
    if t(4,i)==2                                        % PCB
        D = hooke(ptype,Epcb,vpcb);
        Deps0 = Epcb*alphapcb*DeltaT*[1;1:1;0]./(1-2*vpcb);
    end
    if t(4,i)==3                                        % SMD
        D = hooke(ptype,Esmd,vsmd);
        Deps0 = Esmd*alphasmd*DeltaT*[1;1:1;0]./(1-2*vsmd);
    end
    Ke=plante(ex,ey,ep,D);
    fe=plantf(ex,ey,ep,Deps0');
    
    [K,f]=assem(nedof(i,:),K,Ke,f,fe);
end


%% Calculate the boundary Bc.
Bc = zeros(30,1);                                        % Our boundary matrix which will be used in solveq 
counter = 1;
[~, length] = size(e);

for i = 1:length;
    if isMemberL4(i,e,p) || isMemberL6(i,e,p)
        Bc(counter,:) = [2*e(1,i)- 1];               
        counter = counter + 1;
        Bc(counter,:) = [2*e(2,i)- 1];
        counter = counter + 1;
    end
    if isMemberL5(i,e,p)
        Bc(counter,:) = [2*e(1,i)];
        counter = counter + 1;
        Bc(counter,:) = [2*e(2,i)];
        counter = counter + 1;
    end
end
Bc = unique(Bc);
[length, ~] = size(Bc);
Bc = [Bc zeros(length,1)];

%% Solve and plots the displacement in all the nodes.
u = solveq(K,f,Bc);

ed = extract(nedof,u);

%eexpanded = [ed(:,1) ed(:,2) ed(:,3) ed(:,4) ed(:,5) ed(:,6);
%   -ed(:,1) ed(:,2) -ed(:,3) ed(:,4) -ed(:,5) ed(:,6)];    % This is to mirror the solution
%Ex = [Ex;-Ex];                                              % This is to mirror the solution
%Ey = [Ey;Ey];                                               % This is to mirror the solution

plotpar = [1,4,1];
magnfac = 1000;

%eldisp2(Ex,Ey,eexpanded,plotpar,magnfac);                    % Plots the mirrored displacement u
eldisp2(Ex,Ey,ed,plotpar,magnfac);                            % Plots the displacement u
xlabel('width [meter]', 'fontweight','bold','fontsize',15);
ylabel('height [meter]', 'fontweight','bold','fontsize',15);
title('The displacement [Magnification factor 1000]','fontsize',16);
set(gca,'fontsize',12);
%grid on


%% Computes the stresses and strains for each element aswell as the von Mises stress per element

[length,~]=size(nedof);
stress = zeros(nelm,4);
Seff_el = zeros(nelm,1);
for i = 1:length
    ex = Ex(i,:);
    ey = Ey(i,:);
    ep = [ptype 1];
    ed = [u(nedof(i,2)) u(nedof(i,3)) u(nedof(i,4)) u(nedof(i,5)) u(nedof(i,6)) u(nedof(i,7))];     
    % ed: All the displacements in one element
    if t(4,i)==1                            % Solder
        D = hooke(ptype,Esol,vsol);
    end
    if t(4,i)==2                            % PCB 
        D = hooke(ptype,Epcb,vpcb);
    end
    if t(4,i)==3                            % SMD 
        D = hooke(ptype,Esmd,vsmd);
    end
    
    [es,~]=plants(ex,ey,ep,D,ed);
    stress(i,:) = stress(i,:) + es; 
    
    Seff_el(i) = sqrt(stress(i,1)^2 + stress(i,2)^2 + stress(i,3)^2 - stress(i,1)*stress(i,2) - stress(i,1)*stress(i,3) - stress(i,2)*stress(i,3) + 3*stress(i,4)^2);
    % This is the von Mises stress field for an element
end


%% Find the von Mises stress for each node
[length,~] = size(coord);
Seff_nod = zeros(length,1);

for i=1:size(coord,1)
    [c0,c1]=find(edof(:,2:4)==i);
    Seff_nod(i,1)=sum(Seff_el(c0))/size(c0,1);
end

%% 
ed = extract(edof,Seff_nod);
%ed = [ed;ed];                                               % This is to mirror the solution
%Ex = [Ex;-Ex];                                              % This is to mirror the solution
%Ey = [Ey;Ey];                                               % This is to mirror the solution
fill(Ex',Ey',ed','EdgeColor','none')  
xlabel('width [meter]', 'fontweight','bold','fontsize',15);
ylabel('height [meter]', 'fontweight','bold','fontsize',15);
title('The von Mises stress field [N/m^2]','fontsize',16);
set(gca,'fontsize',12);
colorbar;