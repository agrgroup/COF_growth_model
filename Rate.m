function dCdt = Rate(t,C,k)

% This function calculates the rate of the different reactions involved in
% the process of COF growth
% Input: t: time in hours, C: concentration of the species, k: rate
% constants

% Total 17 species are present
% M1 exists in four forms: M1_3, M1_2, M1_1 and M1_0, according to the 
% number of inhibitors attached
M1data=C(1:4);                                                              
M1 = @(x) M1data(4-x);

% M2 can exists as M2_3, M2_2, M2_1 and M2_0, same as M1
M2data=C(5:8);
M2 = @(x) M2data(4-x);

% I1 and I2 are the two monofunctional inhibitors which can react and form
% a compound I12
I = C(9:11);
I1 = I(1);
I2 = I(2);
I12 = I(3);

% Polymers are categorised according to the number of free ends it consists
% and which monomer is attached at the end: P1_2, P1_1, P2_2 and P2_1
P1data = C(12:14);
P1 = @(x) P1data(3-x);
P2data=C(15:17);
P2 = @(x) P2data(3-x);

% km: Rate for monomer attachment, kmr: Rate for monomer de-attachment, ki:
% Rate for inhibitor attachment, kir: Rate for inhibitor de-attachment, Ff:
% Enhancement factor (Forward reaction), Fb: De-enhancement factor
% (Backward reaction)
km = k(1); kmr = k(2); ki = k(3); kir = k(4); Ff = k(5); Fb = k(6);

% Initialisation of the rate matrix
dCdt=zeros(17,1);

%% Equations 1 and 2
dCdt(1) = -ki*M1(3)*I2 + kir*M1(2);                                         % M1_3 + I2 ⇌ M1_2 
dCdt(2) = ki*M1(3)*I2 - kir*M1(2) - ki*M1(2)*I2 + kir*M1(1);                % M1_2 + I2 ⇌ M1_1 and M1_3 + I2 ⇌ M1_2 
dCdt(3) = ki*M1(2)*I2 - kir*M1(1) - ki*M1(1)*I2 + kir*M1(0);                % M1_1 + I2 ⇌ M1_0 and M1_2 + I2 ⇌ M1_1
dCdt(4) = ki*M1(1)*I2 - kir*M1(0);                                          % M1_0 ⇌ M1_1 + I2
dCdt(5) = -ki*M2(3)*I1 + kir*M2(2);                                         % M2_3 + I1 ⇌ M2_2
dCdt(6) = ki*M2(3)*I1 - kir*M2(2) - ki*M2(2)*I1 + kir*M2(1);                % M2_2 + I1 ⇌ M2_1 and M2_3 + I1 ⇌ M2_2
dCdt(7) = ki*M2(2)*I1 - kir*M2(1) - ki*M2(1)*I1 + kir*M2(0);                % M2_1 + I1 ⇌ M2_0 and M2_2 + I1 ⇌ M2_1
dCdt(8) = ki*M2(1)*I1 - kir*M2(0);                                          % M2_0 ⇌ M2_1 + I1
dCdt(9) = -ki*I1*sum(M2(1:3)) + kir*sum(M2(0:2));                           % I1 + M2_3 or M2_2 or M2_1 ⇌ M2_2 or M2_1 or M2_0
dCdt(10) = -ki*I2*sum(M1(1:3)) + kir*sum(M1(0:2));                          % I2 + M1_3 or M1_2 or M1_1 ⇌ M1_2 or M1_1 or M1_0

%% Equation 3

K3 = F3calc();
% The columns for K3 are: 1. Number of free neighbors of M1 (x), 2. Number of 
% free neighbors of M2 (y), 3. Number of P1_1 neighbors (a1), 4.Number of P1_2 
% neighbors (a2), 5. Number of P2_1 neighbors (b1), 6. Number of P2_2
% neighbors (b2)

for i = 1 : length(K3)

 Var = K3(i,:); 
 x = Var(1); y = Var(2); a1 = Var(3); a2 = Var(4); b1 = Var(5); b2 = Var(6);                                  
 
 % Calculating the corresponding forward and backward rates
 Forw3 = km*Ff^(a1+a2+b1+b2)*M1(x)*M2(y)*P1(1)^a1*P1(2)^a2*P2(1)^b1*P2(2)^b2;
 Back3 = kmr*Fb^(-(a1+a2+b1+b2))*P1(x-1-b1-b2)*P2(y-1-a1-a2)*P1(0)^a1*P1(1)^a2*P2(0)^b1*P2(1)^b2;
 
 % Modifying the concentrations 
 dCdt(4-x) = dCdt(4-x) - Forw3 + Back3;                                     % M1_x + M2_y ⇌ P1_(x-1) + P2_(y-1) 
 dCdt(8-y) = dCdt(8-y) - Forw3 + Back3;                                     % M2_y + M1_x ⇌ P2_(y-1) + P1_(x-1)
 dCdt(15-x-b1-b2) = dCdt(15-x-b1-b2) + Forw3 - Back3;                       % P1 
 dCdt(18-y-a1-a2) = dCdt(18-y-a1-a2) + Forw3 - Back3;                       % P2 
 
 % If a2 is greater than 0
 % P1_2 consumed
 dCdt(12) = dCdt(12) + a2*(-Forw3 + Back3);                                 % M1_x + M2_y + a2 * P1_2 ⇌ P1_(x-1) + P2_(y-1) + a2 * P1_1                        
 % P1_1 formed from P1_2
 dCdt(13) = dCdt(13) + a2*(Forw3 - Back3);                                  
 
 % If a1 is greater than 0
 % P1_1 consumed
 dCdt(13) = dCdt(13) + a1*(-Forw3 + Back3);                                 % M1_x + M2_y + a1 * P1_1 ⇌ P1_(x-1) + P2_(y-1) + a1 * P1_0 
 % P1_0 formed from P1_1
 dCdt(14) = dCdt(14) + a1*(Forw3 - Back3);                                  

 % If b2 is greater than 0
 % P2_2 consumed
 dCdt(15) = dCdt(15) + b2*(-Forw3 + Back3);                                 % M1_x + M2_y + b2 * P2_2 ⇌ P1_(x-1) + P2_(y-1) + b2 * P2_1 
 %P2_1 formed from P2_2
 dCdt(16) = dCdt(16) + b2*(Forw3 - Back3);                                  
 
 %b1
 % P2_1 consumed
 dCdt(16) = dCdt(16) + b1*(-Forw3 + Back3);                                 % M1_x + M2_y + b1 * P2_1 ⇌ P1_(x-1) + P2_(y-1) + b1 * P2_0 
 % P2_0 formed from P2_1
 dCdt(17) = dCdt(17) + b1*(Forw3 - Back3);                                 
              
end

%% Equations 4 and 5

K4 = F45calc();
% The columns of K4 are: 1. Number of neighbors of M1 (x), 2. Number of
% neighbors of P2_1 (b1), 3. Number of neighbors of P2_2 (b2)

for i = 1 : length(K4)
  
    Var=K4(i,:); x = Var(1); b1=Var(2); b2 = Var(3);
    
    Forw4 = km*Ff^(b1+b2-1)*M1(x)*P2(1)^b1*P2(2)^b2;
    Back4 = kmr*Fb^(-(b1+b2-1))*P1(x-b1-b2)*P2(0)^b1*P2(1)^b2;
    
    dCdt(4-x) = dCdt(4-x) - Forw4 + Back4;                                  % M1_x ⇌ P1_(x-1)  
 
    dCdt(16-x-b1-b2) = dCdt(16-x-b1-b2) + Forw4 - Back4;                    
 
    % If b2 is greater than 0
    %P2_2 consumed
    dCdt(15) = dCdt(15) + b2*(-Forw4 + Back4);                              % M1_x + b2 * P2_2 ⇌ P1_(x-1) + b2 * P2_1                 
    %P2_1 formed from P2_2
    dCdt(16) = dCdt(16) + b2*(Forw4 - Back4);                              
    
    % If b1 is greater than 0
    %P2_1 consumed
    dCdt(16) = dCdt(16) + b1*(-Forw4 + Back4);                              % M1_x + b1 * P2_1 ⇌ P1_(x-1) + b1 * P2_0
    %P2_0 formed from P2_1
    dCdt(17) = dCdt(17) + b1*(Forw4 - Back4);                               
         
end

K5 = F45calc();
% The columns of K4 are: 1. Number of neighbors of M2, 2. Number of
% neighbors of P1_1, 3. Number of neighbors of P1_2

for i = 1 : length(K5)
    
    Var = K5(i,:); y = Var(1); a1 = Var(2); a2 = Var(3);
    
    Forw5 = km*Ff^(a1+a2-1)*M2(y)*P1(1)^a1*P1(2)^a2;
    Back5 = kmr*Fb^(-(a1+a2-1))*P2(y-a1-a2)*P1(0)^a1*P1(1)^a2;
    
    dCdt(8-y) = dCdt(8-y) - Forw5 + Back5;                                  % M2_y ⇌ P2_y-1
    
    dCdt(19-y-a1-a2) = dCdt(19-y-a1-a2) + Forw5 - Back5;                   
    
    % If a2 is greater than zero
    %P1_2 consumed
    dCdt(12) = dCdt(12) + a2*(-Forw5 + Back5);                              % M2_y + a2 * P1_2 ⇌ P2_y-1 + a2 * P1_1
    %P1_1 formed from P1_2
    dCdt(13) = dCdt(13) + a2*(Forw5 -  Back5);                              
    
    % If a1 is greater than zero
    %P1_1 consumed
    dCdt(13) = dCdt(13) + a1*(-Forw5 + Back5);                              % M2_y + a1 * P1_1 ⇌ P2_y-1 + a1 * P1_0
    %P1_0 formed from P1_1
    dCdt(14) = dCdt(14) + a1*(Forw5 - Back5);
    
end


%% Equations 6 and 7

dCdt(12) = dCdt(12) - ki*P1(2)*I2 + kir*P1(1);                              % P1_2 + I2 ⇌ P1_1
dCdt(13) = dCdt(13) + ki*P1(2)*I2 - kir*P1(1) - ki*P1(1)*I2 + kir*P1(0);    % P1_1 + I2 ⇌ P1_0 and P1_2 + I2 ⇌ P1_1
dCdt(14) = dCdt(14) + ki*P1(1)*I2 - kir*P1(0);                              % P1_0 ⇌ P1_1 + I2

dCdt(10) = dCdt(10) - ki*P1(2)*I2 + kir*P1(1) - ki*P1(1)*I2 + kir*P1(0);    % I2 + P1_x ⇌ P1_(x-1)

dCdt(15) = dCdt(15) - ki*P2(2)*I1 + kir*P2(1);                              % P2_2 + I1 ⇌ P2_1
dCdt(16) = dCdt(16) + ki*P2(2)*I1 - kir*P2(1) - ki*P2(1)*I1 + kir*P2(0);    % P2_1 + I1 ⇌ P2_0 or P2_2 + I1 ⇌ P2_1
dCdt(17) = dCdt(17) + ki*P2(1)*I1 - kir*P2(0);                              % P2_0 ⇌ P2_1 + I1

dCdt(9) = dCdt(9) - ki*P2(2)*I1 + kir*P2(1) - ki*P2(1)*I1 + kir*P2(0);      % I1 + P2_y ⇌ P2_(y-1)

%% Equation 8

dCdt(9) = dCdt(9) - ki*I1*I2 + kir*I12;                                     % I1 + I2 ⇌ I1I2
dCdt(10) = dCdt(10) - ki*I2*I1 + kir*I12;                                   % I2 + I1 ⇌ I1I2
dCdt(11) = dCdt(11) + ki*I1*I2 - kir*I12;                                   % I1I2 ⇌ I1 + I2
