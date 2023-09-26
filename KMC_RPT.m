function [tf,tknock_timeseries,sites_list,C] = KMC_RPT(basedir)

% This function performs a single kinetic Monte Carlo simulation for the
% COF growth in a 2D hexagonal lattice 

global basedirectory N Lx Ly Lz NA area_site lattice_const ...
    num_M1_3 num_M2_3 num_I1 num_I2 km ki kir F_f P;

% Initialize time
t0 = 0;

% Set base directory global variable
basedirectory = basedir;

% Set grid size
N = 50;

% Set lattice size
lattice_const = 18.72115;
area_site = 0.50*sqrt(3)*(lattice_const^2)*1e-20;

box_x = N*2*(lattice_const/sqrt(3))*cos(pi/6);
box_y = N*(lattice_const/sqrt(3))*(2+2*cos(pi/3))/2;

% Lz is calculated by matching the volume used in the experiments for COF growth 
Lx = box_x; Ly = box_y; Lz = 1094;


% Values of physical constant
NA = 6.022e23;    % Avogadro's number in mol^-1

% Initial Concentrations of monomers and inhibitors
M1_3 = 0.005;  %mmol/mL
M2_3 = 0.005;
I_1 = 12*M1_3; %12C Case
I_2 = 12*M2_3;

% Number of Molecules calculation

num_M1_3 = Lx*Ly*Lz*M1_3*1e3*1e-30*NA;
num_M2_3 = Lx*Ly*Lz*M2_3*1e3*1e-30*NA;
num_I1 = Lx*Ly*Lz*I_1*1e3*1e-30*NA;
num_I2 = Lx*Ly*Lz*I_2*1e3*1e-30*NA;

P = 1; %activity

% Rate Constants (hr^-1) obtained by theoretical model
km = 302.1246*20;
kir = 4.6879*0.03;
ki = 1301.8;

% Enhancement Factor for forward reaction
F_f = 1.3585;

% Initialize 2D material layer without any P atoms
C = zeros(N,N,2);
num_atoms = [0,0];

% Do the KMC algorithm
tic
[tf,tknock_timeseries,sites_list,C] = kmc(C,num_atoms,t0,N);
toc

end

function [t,tknock_timeseries,sites_list,C] = kmc(C0,num_atoms,t0,N)
% This function performs the KMC algorithm and returns final configuration
% of the system

% Inputs: initial surface configuration (C0), number of atoms on COF
% sublattices (num_atoms), initial time (t0) and grid size (N)
global num_M1_3 num_M2_3 num_I1 num_I2 Lx Ly Lz NA F_f P km ki;

t = t0; % initial time
count = 0; % counter for number of time steps
C = C0; % store initial state of system in variable C

% Calculate initial rates
n_Species = [num_M1_3 0 0 0 num_M2_3 0 0 0 num_I1 num_I2 0];

% Initialise pairs list for set 2 reactions
pairs_list = pairs(C,N,n_Species);

% Assuming one of the set 2 reactions was picked
% Initialize the sites list and choose a pair for P1, P2 addition randomly
sites_list = [];

num_nucleate = 1;
for nucleate = 1 : num_nucleate
    list=[];
    for i = 1:length(pairs_list)
        
        % There are numerous possible pairs for P1-P2 addition for a particular
        % reaction
        if (pairs_list(i,5)==15 && pairs_list(i,6)>0)
            list = [list i];
        end
    end

% Pick a particular pair randomly
sel = randi(length(list));
pair_sel = pairs_list(list(sel),1:4);
reactype_sel = pairs_list(list(sel),5);
P1P2_type = Ptype(15);

% Add the P sites to the sites_list and update the pairs_list
[sites_list,C,num_atoms,occupancy] = sites(C,P1P2_type,num_atoms,sites_list,pair_sel,n_Species);
[pairs_list,rows_updated] = pairslist_update_reac2(pairs_list,pair_sel,occupancy);

% Upating the n_Species matrix
n_Species = nChange(reactype_sel,n_Species);

end

num_dt=sum(num_atoms);
tknock_timeseries = zeros(num_dt+1,1);

% Combining the rates arrays for three set of reactions
Rates_1 = reactions(n_Species); %Rates for 1st set reactions
Rates_2 = pairs_list(:,6); %Rates for 2nd set reactions
Rates_3 = sites_list(:,5); %Rates for 3rd set reactions

Rates = [Rates_1(1:14); Rates_2; Rates_3]; %Concatenate all the rates in one array

ReactionSel = []; %To save the reaction number picked at each iteration
Time_total = [];
Monomer1_3 = []; Monomer1_2 = []; Monomer1_1 = []; Monomer1_0 = [];
Monomer2_3 = []; Monomer2_2 = []; Monomer2_1 = []; Monomer2_0 = [];
Inhibitor1 = [];
Inhibitor2 = [];
Yld = [];

%Main Loop
z = 0;

while t < 24
    z = z + 1;

    % Generate two uniformly-distributed random numbers u1 and u2
    random = rand(1,2);
    u1=random(1);
    u2=random(2);
    
    % Compute effective rate vector
    [num,~] = size(Rates);
    rcum = cumsum(Rates);
    rcum = rcum/rcum(num); % normalize with total rate
    
    tknock_timeseries(count+1) = 1/sum(Rates);
    % Find which of the sites are going to have a reaction, based on
    % Gillespie's algorithm
    % disp(sites_list);
    reac_sel=0;
    for i=1:num
        if (i==1 && u1<=rcum(i))
            reac_sel=i;
        elseif (i>1 && u1>rcum(i-1) && u1<=rcum(i))
            reac_sel=i;
        end
    end
    
    % Increment the system time based on a Poisson process
    rnet = sum(Rates);
    if (rnet > 0)
        dt = log(1/u2)/rnet;
       % dt = 1/rnet;
        t = t + dt;
        count = count + 1;
    end
    
    % 1st set of reaction: update the concentrations 
    if (reac_sel <= 14)
        
        n_Species = nChange(reac_sel,n_Species);
        disp([z reac_sel]);
        
        % 2nd set of reaction: Add the pair to the sites_list and update the pairs_list
    elseif (reac_sel >= 15 && reac_sel <= length(pairs_list)+14)
        
        pair_sel = pairs_list(reac_sel-14,1:4); %y1, x1, y2, x2
        reactype_sel = pairs_list(reac_sel-14,5); %Reaction number
        P1P2_type = Ptype(reactype_sel); %What kind of P1 and P2 to add on the lattice sites
        
        %Adding the sites to the sites_list and updating the pairs_list
        [sites_list,C,num_atoms] = sites(C,P1P2_type,num_atoms,sites_list,pair_sel,n_Species);
        [pairs_list,~] = pairslist_update_reac2(pairs_list,pair_sel,occupancy);
        
        %Updating the n_Species
        n_Species = nChange(reactype_sel,n_Species);
        
        reac_sel = reactype_sel;
        disp([z reac_sel]);
        
        % 3rd set of reactions: Update the sites_list and the pairs_list, all the rates
    elseif (reac_sel > length(pairs_list) + 14)
        
        modified_sites = [];
        
        % Find (r,c) position of the lattice site selected to be formed
        sel_row = reac_sel - length(pairs_list) - 14;
        pos_sel_r = sites_list(sel_row,2);
        pos_sel_c = sites_list(sel_row,3);
        
        % Find sublattice at which atom is going to attach
        sublattice_sel = sites_list(sel_row,1);
        modified_sites = [modified_sites; pos_sel_r pos_sel_c sublattice_sel];
        
        % Find if the selected lattice site is already occupied
        site_status = sites_list(sel_row,4);
        inhibitor_status = sites_list(sel_row,6);
        
        %Find what kind of reaction is to be picked
        if (site_status == 1 && inhibitor_status == 3) % Inhibitor to be removed
            if (sublattice_sel == 1)
                reac_type = 33;
            else
                reac_type = 31;
            end
            
            C(pos_sel_r,pos_sel_c,sublattice_sel) = 0;
            num_atoms(sublattice_sel) = num_atoms(sublattice_sel) - 1;
            reac_sel = reac_type;
            
        elseif (site_status == 0 && inhibitor_status == 0) % Monomer or inhibitor to be added
            
            Species = n_Species./(Lx*Ly*Lz*1e3*1e-30*NA);
            % Species = [M1_3 M1_2 M1_1 M1_0 M2_3 M2_2 M2_1 M2_0 I_1 I_2 I1I2];
            M1data=Species(1:4);
            M1 = @(x) M1data(4-x);
            M2data=Species(5:8);
            M2 = @(x) M2data(4-x);
            I=Species(9:11);
            
            % Storing all the possible rates
            sites_rates1 = [km*F_f^2*M1(3) km*F_f^2*M1(2) km*F_f^1*P*M1(1) ki*I(1)];
            sites_rates2 = [km*F_f^2*M2(3) km*F_f^2*M2(2) km*F_f^1*P*M2(1) ki*I(2)];
            sites_rates = [sites_rates1; sites_rates2];
            
            reac_type1 = [24 25 26 32];
            reac_type2 = [27 28 29 30];
            reac_type = [reac_type1; reac_type2];
            
            num_neighbor = count_neighbors(pos_sel_r, pos_sel_c, sublattice_sel,C);
            
            % Checking if addition of M(1) or M(2) is possible (as they have inhibitors), if near_neighbors are not occupied
            if num_neighbor == 3
                sites_rates(sublattice_sel,2:3)=0;
            elseif num_neighbor == 2
                sites_rates(sublattice_sel,3)=0;
            end
            
            %Performing Gillespie's algorithm to pick a reaction for the
            %particular site
            v=random(1);
            num1=length(sites_rates(sublattice_sel,:));
            rcum1 = cumsum(sites_rates(sublattice_sel,:));
            rcum1=rcum1/rcum1(num1); % normalize with total rate
            
            reac_selected = 0;
            for i = 1 : num1
                if (i == 1 && v <= rcum1(i))
                    reac_selected = i;
                elseif (i > 1 && v > rcum1(i-1) && v <= rcum1(i))
                    reac_selected = i;
                end
                
            end
            reac_type = reac_type(sublattice_sel,reac_selected);
            inhibitor = reac_selected - 1; %Number of inhibitors to be added at the lattice sites (for M(3) = 0, M(2) = 1, M(1) = 2, for inhibitor nothing)
            
            if inhibitor == 0 %M(3) addition
                % If selected lattice site is 'vacant', then update site
                % status to 'occupied'
                C(pos_sel_r,pos_sel_c,sublattice_sel) = sublattice_sel;
                num_atoms(sublattice_sel) = num_atoms(sublattice_sel) + 1;
                
            elseif inhibitor == 3 %I addition
                % If selected lattice site is 'vacant', then update site
                % status to 'occupied' and 'inhibitor'
                C(pos_sel_r,pos_sel_c,sublattice_sel) = sublattice_sel + 2;
                num_atoms(sublattice_sel) = num_atoms(sublattice_sel) + 1;
                
            elseif inhibitor == 1 %M(2) addition
                C(pos_sel_r,pos_sel_c,sublattice_sel) = sublattice_sel;
                num_atoms(sublattice_sel) = num_atoms(sublattice_sel) + 1;
                
                % Inhibitor addition in the neighbor sites
                [nearest,~] = neighbor_sites(pos_sel_r,pos_sel_c,sublattice_sel);
                
                % Adding one inhibitor
                countinhib = 0;
                for i = 1 : size(nearest)
                    if (C(nearest(i,1),nearest(i,2),nearest(i,3)) == 0)
                        C(nearest(i,1),nearest(i,2),nearest(i,3)) = nearest(i,3) + 2;
                        num_atoms(nearest(i,3)) = num_atoms(nearest(i,3))+1;
                        countinhib = countinhib + 1;
                        modified_sites = [modified_sites; nearest(i,1),nearest(i,2),nearest(i,3)];
                        if (countinhib==1)
                            break;
                        end
                    end
                end
                
            elseif inhibitor == 2 %M(1) addition
                C(pos_sel_r,pos_sel_c,sublattice_sel) = sublattice_sel;
                num_atoms(sublattice_sel) = num_atoms(sublattice_sel)+1;
                
                %Inhibitors addition in the neighbor sites
                [nearest,~]=neighbor_sites(pos_sel_r,pos_sel_c,sublattice_sel);
                
                % Additing two inhibitors
                countinhib = 0;
                for i = 1 : size(nearest)
                    if (C(nearest(i,1),nearest(i,2),nearest(i,3)) == 0)
                        C(nearest(i,1),nearest(i,2),nearest(i,3)) = nearest(i,3)+2;
                        num_atoms(nearest(i,3)) = num_atoms(nearest(i,3))+1;
                        countinhib = countinhib + 1;
                        modified_sites = [modified_sites; nearest(i,1),nearest(i,2),nearest(i,3)];
                    end
                end
                
            end
            
        end

        % Update sites_list
        [sites_list, ~] = sites_change(C,sublattice_sel,modified_sites,sites_list,site_status,inhibitor,n_Species);        
        
        % Update pairs_list
        pairs_list = pairslist_update_reac3(pairs_list,modified_sites,site_status);
        
        % Which reaction number was picked
        reac_sel = reac_type;
        
        % Update the n_Species array and the rates in pairs_list and
        % sites_list accordingly
        n_Species = nChange(reac_sel,n_Species);
        disp([z reac_sel]);        
      
    end

    % Calculation of yield
    Species = n_Species./(Lx*Ly*Lz*1e3*1e-30*NA);
    
    Yield = (0.005 - Species(1) - Species(2) - Species(3) - Species(4))/0.005;
    disp(t); disp(Yield);
   
    Rates_1 = reactions(n_Species);
    Rates_2 = pairs_list(:,6);
    Rates_3 = sites_list(:,5);
    
    Rates = [Rates_1(1:14); Rates_2; Rates_3];
    
    % Visualization
    if (mod(z,250)==0)
        Species = n_Species./(Lx*Ly*Lz*1e3*1e-30*NA);
        visualize(C,N,z,reac_sel,Species);
        
    end     
    %Storing reaction number and concentrations after each iteration

    ReactionSel = [ReactionSel reac_sel];
    Time_total =[Time_total t];
    Monomer1_3 = [Monomer1_3 Species(1)];
    Monomer1_2 = [Monomer1_2 Species(2)];
    Monomer1_1 = [Monomer1_1 Species(3)];
    Monomer1_0 = [Monomer1_0 Species(4)];
    Monomer2_3 = [Monomer2_3 Species(5)];
    Monomer2_2 = [Monomer2_2 Species(6)];
    Monomer2_1 = [Monomer2_1 Species(7)];
    Monomer2_0 = [Monomer2_0 Species(8)];
    Inhibitor1 = [Inhibitor1 Species(9)];
    Inhibitor2 = [Inhibitor2 Species(10)];

    Yld = [Yld Yield];

end

save('ReactionSel_18C_final.mat','ReactionSel','Time_total','Monomer1_3','Monomer1_2',...
'Monomer1_1','Monomer1_0','Monomer2_3','Monomer2_2','Monomer2_1','Monomer2_0','Inhibitor1','Inhibitor2','Yld');
end

function nNew = nChange(reac_sel,n_Species)

%This function changes the number of molecules according to the reaction
%picked in the n_Species array

%Input: reac_sel: Reaction number picked, n_Species: Array of number of
%molecules of each species

num_M1 = n_Species(1:4); num_M2 = n_Species(5:8); nI = n_Species(9:11);
nM1 = @(x) num_M1(4-x);
nM2 = @(x) num_M2(4-x);

n_Change = [nM1(3)-1 nM1(2)+1 nM1(1) nM1(0) nM2(3) nM2(2) nM2(1) nM2(0) nI(1) nI(2)-1 nI(3); % M1(3) + I(2) -> M1(2) %% Set1
    nM1(3)+1 nM1(2)-1 nM1(1) nM1(0) nM2(3) nM2(2) nM2(1) nM2(0) nI(1) nI(2)+1 nI(3);         % M1(2) -> M1(3) +I(2)
    nM1(3) nM1(2)-1 nM1(1)+1 nM1(0) nM2(3) nM2(2) nM2(1) nM2(0) nI(1) nI(2)-1 nI(3);         % M1(2) + I(2) -> M1(1)
    nM1(3) nM1(2)+1 nM1(1)-1 nM1(0) nM2(3) nM2(2) nM2(1) nM2(0) nI(1) nI(2)+1 nI(3);         % M1(1) -> M1(2) +I(2)
    nM1(3) nM1(2) nM1(1)-1 nM1(0)+1 nM2(3) nM2(2) nM2(1) nM2(0) nI(1) nI(2)-1 nI(3);         % M1(1) + I(2) -> M1(0)
    nM1(3) nM1(2) nM1(1)+1 nM1(0)-1 nM2(3) nM2(2) nM2(1) nM2(0) nI(1) nI(2)+1 nI(3);         % M1(0) -> M1(1) + I(2)
    nM1(3) nM1(2) nM1(1) nM1(0) nM2(3)-1 nM2(2)+1 nM2(1) nM2(0) nI(1)-1 nI(2) nI(3);         % M2(3) + I(1) -> M2(2)
    nM1(3) nM1(2) nM1(1) nM1(0) nM2(3)+1 nM2(2)-1 nM2(1) nM2(0) nI(1)+1 nI(2) nI(3);         % M2(2) -> M2(3) + I(1)
    nM1(3) nM1(2) nM1(1) nM1(0) nM2(3) nM2(2)-1 nM2(1)+1 nM2(0) nI(1)-1 nI(2) nI(3);         % M2(2) + I(1) -> M2(1)
    nM1(3) nM1(2) nM1(1) nM1(0) nM2(3) nM2(2)+1 nM2(1)-1 nM2(0) nI(1)+1 nI(2) nI(3);         % M2(1) -> M2(2) + I(1)
    nM1(3) nM1(2) nM1(1) nM1(0) nM2(3) nM2(2) nM2(1)-1 nM2(0)+1 nI(1)-1 nI(2) nI(3);         % M1(1) + I(1) -> M1(0)
    nM1(3) nM1(2) nM1(1) nM1(0) nM2(3) nM2(2) nM2(1)+1 nM2(0)-1 nI(1)+1 nI(2) nI(3);         % M1(0) -> M1(1) + I(1)
    nM1(3) nM1(2) nM1(1) nM1(0) nM2(3) nM2(2) nM2(1) nM2(0) nI(1)-1 nI(2)-1 nI(3)+1;         % I(1) + I(2) -> I1I2
    nM1(3) nM1(2) nM1(1) nM1(0) nM2(3) nM2(2) nM2(1) nM2(0) nI(1)+1 nI(2)+1 nI(3)-1;         % I1I2 -> I(1) + I(2)
    nM1(3)-1 nM1(2) nM1(1) nM1(0) nM2(3)-1 nM2(2) nM2(1) nM2(0) nI(1) nI(2) nI(3);           % M1(3) + M2(3) +.. -> P1(2) + P2(2) +.. %%Set2
    nM1(3) nM1(2)-1 nM1(1) nM1(0) nM2(3)-1 nM2(2) nM2(1) nM2(0) nI(1) nI(2) nI(3);           % M1(2) + M2(3) +.. -> P1(1) + P2(2) +..
    nM1(3) nM1(2) nM1(1)-1 nM1(0) nM2(3)-1 nM2(2) nM2(1) nM2(0) nI(1) nI(2) nI(3);           % M1(1) + M2(3) +.. -> P1(0) + P2(2) +..
    nM1(3)-1 nM1(2) nM1(1) nM1(0) nM2(3) nM2(2)-1 nM2(1) nM2(0) nI(1) nI(2) nI(3);           % M1(3) + M2(2) +.. -> P1(2) + P2(1) +..
    nM1(3) nM1(2)-1 nM1(1) nM1(0) nM2(3) nM2(2)-1 nM2(1) nM2(0) nI(1) nI(2) nI(3);           % M1(2) + M2(2) +.. -> P1(1) + P2(1) +..
    nM1(3) nM1(2) nM1(1)-1 nM1(0) nM2(3) nM2(2)-1 nM2(1) nM2(0) nI(1) nI(2) nI(3);           % M1(1) + M2(2) +.. -> P1(0) + P2(1) +..
    nM1(3)-1 nM1(2) nM1(1) nM1(0) nM2(3) nM2(2) nM2(1)-1 nM2(0) nI(1) nI(2) nI(3);           % M1(3) + M2(1) +.. -> P1(2) + P2(0) +..
    nM1(3) nM1(2)-1 nM1(1) nM1(0) nM2(3) nM2(2) nM2(1)-1 nM2(0) nI(1) nI(2) nI(3);           % M1(2) + M2(1) +.. -> P1(1) + P2(0) +..
    nM1(3) nM1(2) nM1(1)-1 nM1(0) nM2(3) nM2(2) nM2(1)-1 nM2(0) nI(1) nI(2) nI(3);           % M1(1) + M2(1) +.. -> P1(0) + P2(0) +..
    nM1(3)-1 nM1(2) nM1(1) nM1(0) nM2(3) nM2(2) nM2(1) nM2(0) nI(1) nI(2) nI(3);             % M1(3) + .. -> P1(2) +..  %%Set3
    nM1(3) nM1(2)-1 nM1(1) nM1(0) nM2(3) nM2(2) nM2(1) nM2(0) nI(1) nI(2) nI(3);             % M1(2) + .. -> P1(1) +..
    nM1(3) nM1(2) nM1(1)-1 nM1(0) nM2(3) nM2(2) nM2(1) nM2(0) nI(1) nI(2) nI(3);             % M1(1) + .. -> P1(0) +..
    nM1(3) nM1(2) nM1(1) nM1(0) nM2(3)-1 nM2(2) nM2(1) nM2(0) nI(1) nI(2) nI(3);             % M2(3) + .. -> P2(2) +..
    nM1(3) nM1(2) nM1(1) nM1(0) nM2(3) nM2(2)-1 nM2(1) nM2(0) nI(1) nI(2) nI(3);             % M2(2) + .. -> P2(1) +..
    nM1(3) nM1(2) nM1(1) nM1(0) nM2(3) nM2(2) nM2(1)-1 nM2(0) nI(1) nI(2) nI(3);             % M2(1) + .. -> P2(0) +..
    nM1(3) nM1(2) nM1(1) nM1(0) nM2(3) nM2(2) nM2(1) nM2(0) nI(1) nI(2)-1 nI(3);             % P1(x) + I(2) -> P1(x-1)
    nM1(3) nM1(2) nM1(1) nM1(0) nM2(3) nM2(2) nM2(1) nM2(0) nI(1) nI(2)+1 nI(3);             % P1(x-1) -> P1(x) + I(2)
    nM1(3) nM1(2) nM1(1) nM1(0) nM2(3) nM2(2) nM2(1) nM2(0) nI(1)-1 nI(2) nI(3);             % P2(y) + I(1) -> P2(y-1)
    nM1(3) nM1(2) nM1(1) nM1(0) nM2(3) nM2(2) nM2(1) nM2(0) nI(1)+1 nI(2) nI(3)];            % P2(y-1) -> P2(y) + I(1)

nNew = n_Change(reac_sel,:);
end

function Rates = reactions(n_Species)

% This function calculates the Reaction rates
global km ki kir Lx Ly Lz NA P;

Species = n_Species./(Lx*Ly*Lz*1e3*1e-30*NA);

% Species = [M1_3 M1_2 M1_1 M1_0 M2_3 M2_2 M2_1 M2_0 I_1 I_2 I1I2];
M1data = Species(1:4);
M1 = @(x) M1data(4-x);
M2data=Species(5:8);
M2 = @(x) M2data(4-x);
I=Species(9:11);
I1=I(1);
I2=I(2);
I12=I(3);
P2 = P; P1 = P; P0 = P;

Rates = [ki*M1(3)*I2; kir*M1(2); ki*M1(2)*I2; kir*M1(1); ki*M1(1)*I2; kir*M1(0); % Set1
    ki*M2(3)*I1; kir*M2(2); ki*M2(2)*I1; kir*M2(1); ki*M2(1)*I1; kir*M2(0);
    ki*I1*I2; kir*I12;
    km*M1(3)*M2(3); km*M1(2)*M2(3); km*M1(1)*M2(3); % Set 2
    km*M1(3)*M2(2); km*M1(2)*M2(2); km*M1(1)*M2(2);
    km*M1(3)*M2(1); km*M1(2)*M2(1); km*M1(1)*M2(1)];

end

function pairs_list = pairs(C,N,n_Species)

% Function to maintain the reaction rates calculation for set 2 reactions
% Pairs List Columns: 1. y1 2. x1 3. y2 4. x2 5. Reaction no. 6. Rate

global Lx Ly Lz NA km;

Species = n_Species./(Lx*Ly*Lz*1e3*1e-30*NA);

% Species are [M1_3 M1_2 M1_1 M1_0 M2_3 M2_2 M2_1 M2_0 I_1 I_2 I1I2];
M1data=Species(1:4);
M1 = @(x) M1data(4-x);
M2data=Species(5:8);
M2 = @(x) M2data(4-x);

% Calculating rates assuming no P sites are available on the lattices
Rates = [M1(3)*M2(3) M1(2)*M2(3) M1(1)*M2(3) M1(3)*M2(2) M1(2)*M2(2)...
    M1(1)*M2(2) M1(3)*M2(1) M1(2)*M2(1) M1(1)*M2(1)];
Rates = Rates*km;

pairs_list=[];
for i=1:N
    for j=1:N
        
        % To avoid repition of pairs, loop is run according to only one
        % sublattice (1)
        [nearest,~] = neighbor_sites(i,j,1);
        
        % Each row is added nine times, as nine possible reactions for each
        % pair
        p1 = [i j nearest(1,1) nearest(1,2)];
        pairs_list = [pairs_list; repmat(p1,9,1)];
        
        p2 = [i j nearest(2,1) nearest(2,2)];
        pairs_list = [pairs_list;
            repmat(p2,9,1)];
        
        p3 = [i j nearest(3,1) nearest(3,2)];
        pairs_list = [pairs_list;
            repmat(p3,9,1)];
        
    end
end

for i = 1:9:length(pairs_list)-1
    
    pairs_list(i:i+8,5)=15:23; %Assigning reaction numbers
    pairs_list(i:i+8,6)=Rates; %Assinging corresponding rates

end

end

function [pairs_list,pairs_list_update] = pairslist_update_reac2(pairs_list,pair_sel,occupancy)

% This function updates the pair list according to the 2nd set reaction picked

global F_f P;
pairs_list_update = [];

% Storing the neighbors of the 2 sites of the selected pair
[near1,~] = neighbor_sites(pair_sel(1),pair_sel(2),1);
[near2,~] = neighbor_sites(pair_sel(3),pair_sel(4),2);

% All the pairs where the selected lattices are present, putting Rate  = 0;
for i = 1 : length(pairs_list)
    if ((pairs_list(i,1)==pair_sel(1) && pairs_list(i,2)==pair_sel(2)) ||...
            (pairs_list(i,3)==pair_sel(3) && pairs_list(i,4)==pair_sel(4)))
        pairs_list(i,6)=0;
        pairs_list_update = [pairs_list_update; pairs_list(i,:)];
    end
end

% Updating the pairs where selected lattice sites are their neighbors
for j = 1: length(near1)
    for i = 1 : length(pairs_list)
        
        % Introducing the enhancement factor
        if ((near1(j,1)==pairs_list(i,3) && near1(j,2)==pairs_list(i,4) && pairs_list(i,6) > 0) ||...
                (near2(j,1)==pairs_list(i,1) && near2(j,2)==pairs_list(i,2) && pairs_list(i,6) > 0))
            
        end
        
        % For reaction involving M1(0) or M2(0), even if 1 neighbor is
        % occupied, this reaction cannot take place
        if ((near1(j,1)==pairs_list(i,3) && near1(j,2)==pairs_list(i,4)) && pairs_list(i,6)>0 && pairs_list(i,5)>20)
            
            pairs_list(i,6)=0;
            
        end
        
        if ((near2(j,1)==pairs_list(i,1) && near2(j,2)==pairs_list(i,2)) && pairs_list(i,6)>0 && (pairs_list(i,5) == 17 || pairs_list(i,5) == 20 || pairs_list(i,5) == 23))
            pairs_list(i,6)=0;
            
        end
    end
end

% If inhibitors were also added to the lattice sites, update the pairs
% rates involving these sites and their corresponding neighbor's also
[occupancyLength,~] = size(occupancy);

if occupancyLength > 2
    for k = 3:occupancyLength
        [near,~]=neighbor_sites(occupancy(k,2),occupancy(k,3),occupancy(k,1));
        for j = 1:length(near)
            for i = 1:length(pairs_list)
                % pairs involving inhibitor sites
                if ((occupancy(k,2)==pairs_list(i,1) && occupancy(k,3) == pairs_list(i,2)) ||...
                        (occupancy(k,2)==pairs_list(i,3) && occupancy(k,3) == pairs_list(i,4)))
                    pairs_list(i,6)=0;
                    pairs_list_update = [pairs_list_update; pairs_list(i,:)];
                end
                
                % For reaction involving M1(0) or M2(0), even if 1 neighbor is
                % occupied, this reaction cannot take place
                if (near(j,1)==pairs_list(i,3) && near(j,2)==pairs_list(i,4) && pairs_list(i,5)>20 ||...
                        near(j,1)==pairs_list(i,1) && near(j,2)==pairs_list(i,2) && pairs_list(i,5)>20)
                    pairs_list(i,6)=0;
                    pairs_list_update = [pairs_list_update; pairs_list(i,:)];
                end
            end
        end
    end
end

end

function pairs_list = pairslist_update_reac3(pairs_list,modified_sites,site_status)

% This function updates the pairs_list whenever 3rd set reaction is picked
% (very similar to the function pairslist_update_reac3)

global F_f P

pair_sel = modified_sites(1,:);

[near,~] = neighbor_sites(pair_sel(1),pair_sel(2),pair_sel(3));

% Used for checking the neighbor's loop
if pair_sel(3)== 1
    k1 = 3; k2 = 4;
else
    k1 = 1; k2 = 2;
end

if site_status == 0
    % All the pairs where the selected lattices are present, putting Rate  = 0;
    for i = 1 : length(pairs_list)
        if ((pairs_list(i,1)==pair_sel(1) && pairs_list(i,2)==pair_sel(2)) || ...
                (pairs_list(i,3)==pair_sel(1) && pairs_list(i,4)==pair_sel(2)))
            pairs_list(i,6)=0;
            
        end
    end
    
    % Updating the pairs where selected lattice sites are their neighbors
    for j = 1: length(near)
        for i = 1 : length(pairs_list)
            %Introducing the enhancement factor
            if (near(j,1) == pairs_list(i,k1) && near(j,2) == pairs_list(i,k2) && pairs_list(i,6)>0)
                pairs_list(i,8)=pairs_list(i,6)*F_f*P;
                
            end
            
            % For reaction involving M1(0) or M2(0), even if 1 neighbor is
            % occupied, this reaction cannot take place
            if pair_sel(3) == 2
                if (near(j,1)==pairs_list(i,k1) && near(j,2)==pairs_list(i,k2) && pairs_list(i,5) > 20)
                    pairs_list(i,6)=0;
                    
                end
                
            else
                if (near(j,1)==pairs_list(i,k1) && near(j,2)==pairs_list(i,k2) &&...
                        (pairs_list(i,5) == 17 || pairs_list(i,5) == 20 || pairs_list(i,5) == 23))
                    pairs_list(i,6)=0;
                end
                
            end
        end
    end
    
    % If inhibitors were also added to the lattice sites, update the pairs
    % rates involving these sites and their corresponding neighbor's also
    [modifiedLength,~] = size(modified_sites);
    
    if modifiedLength > 1
        for k = 2 : modifiedLength
            [near,~] = neighbor_sites(modified_sites(k,1),modified_sites(k,2),modified_sites(k,3));
            
            for j = 1:length(near)
                for i = 1:length(pairs_list)
                    
                    % pairs involving inhibitor sites
                    if (modified_sites(k,1) == pairs_list(i,1) && modified_sites(k,2) == pairs_list(i,2)) ||...
                            (modified_sites(k,1) == pairs_list(i,3) && modified_sites(k,2) == pairs_list(i,4))
                        pairs_list(i,6)=0;
                        
                    end
                    
                    % For reaction involving M1(0) or M2(0), even if 1 neighbor is
                    % occupied, this reaction cannot take place
                    if modified_sites(k,3)==2
                        if (near(j,1)==pairs_list(i,3) && near(j,2)==pairs_list(i,4) && pairs_list(i,5) > 20) ||...
                                (near(j,1)==pairs_list(i,1) && near(j,2)==pairs_list(i,2) && pairs_list(i,5) > 20)
                            pairs_list(i,6)=0;
                            
                        end
                    else
                        if (near(j,1)==pairs_list(i,3) && near(j,2)==pairs_list(i,4) && (pairs_list(i,5) == 17 || pairs_list(i,5) == 20 || pairs_list(i,5) == 23)) ||...
                                (near(j,1)==pairs_list(i,1) && near(j,2)==pairs_list(i,2) && (pairs_list(i,5) == 17 || pairs_list(i,5) == 20 || pairs_list(i,5) == 23))
                            pairs_list(i,6)=0;
                        end
                    end
                end
            end
        end
        
        
    end
    
end
end

function type = Ptype(reac_sel)

%Function gives the number of vacant neighbors of P1 and P2 to know what
%all species to add on the lattice sites in the init function

P1P2_type = [2, 2; 1, 2; 0, 2;
    2, 1; 1, 1; 0, 1;
    2, 0; 1, 0; 0, 0];
type = P1P2_type(reac_sel - 14, :);
end

function [sites_list,C,num_atoms,occupancy]=sites(C,P1P2_type,num_atoms,sites_list,pairs_list,n_Species)

% Compute the positions of all possible relavent interactions
% The columns are 1. sublattice,  2.row-position (y), 3. column-position (x), 4. occupancy flag (1-occupied by atom,
% 0-occupied by vacancy), 5. Rates 6. I/M indicator (1-Monomer, 3-Inhibitor)

global P F_f km ki kir Lx Ly Lz NA;

Species=n_Species./(Lx*Ly*Lz*1e3*1e-30*NA);
% Species = [M1_3 M1_2 M1_1 M1_0 M2_3 M2_2 M2_1 M2_0 I_1 I_2 I1I2];
M1data=Species(1:4);
M1 = @(x) M1data(4-x);
M2data=Species(5:8);
M2 = @(x) M2data(4-x);
I=Species(9:11);
I = [I(2) I(1)];

% Add the P1 P2 atoms on lattice sites with I atoms (if any)
[occupancy,C,num_atoms] = init(C,P1P2_type,num_atoms,pairs_list);
[num,~] = size(occupancy);

[sitesLength,~] = size(sites_list);
[occupancyLength,~] = size(occupancy);

if sitesLength == 0
        
        sites_list(:,[1,2,3]) = occupancy(:,[1,2,3]); % store sublattice, row location, column location
        sites_list(:,4) = 1;
        sites_list(:,5) = 0;
        
        if occupancyLength > 2
            sites_list(3:end,5) = kir*P ;
            
        end
        
        for i = 1 : occupancyLength
            if (C(sites_list(i,2),sites_list(i,3),sites_list(i,1)) == sites_list(i,1))
                sites_list(i,6) = 1;
            else
                sites_list(i,6) = 3;
            end
        end
        modified = 1:1:occupancyLength;
    
end

% Adding the new P sites to the sites_list

if sitesLength > 0
    
modified = [];
for i = 1  : occupancyLength
    [lia,locb] =   ismember([occupancy(i,1), occupancy(i,2), occupancy(i,3)],sites_list(:,[1,2,3]),'rows');
    
    
    if any (~lia)
        indx = sitesLength + 1;
        sites_list(indx,[1,2,3]) = occupancy(i,[1,2,3]); % store sublattice, row location, column location
        sites_list(indx,4) = 1;
        sites_list(indx,5) = 0;
        [sitesLength,~] = size(sites_list);
        
        if i > 2
            sites_list(indx,5) = kir*P ;
            
        end
        
        if (C(sites_list(indx,2),sites_list(indx,3),sites_list(indx,1)) == sites_list(indx,1))
            sites_list(indx,6) = 1;
        else
            sites_list(indx,6) = 3;
        end
        modified = [modified, indx];
        
    end
    
    if locb > 0
        
        sites_list(locb,4) = 1;
        sites_list(locb,5) = 0;
        
        if i > 2
            sites_list(locb,5) = kir*P ;
            
        end
        
        if (C(sites_list(locb,2),sites_list(locb,3),sites_list(locb,1)) == sites_list(locb,1))
            sites_list(locb,6) = 1;
        else
            sites_list(locb,6) = 3;
        end
        
        modified = [modified, locb];
        
    end
    
    
end
end


modifiedLength = length(modified);


% To update number of nearest neighbors and neighbor sites data
for i = 1 : modifiedLength
    
    sublattice = sites_list(modified(i),1);
    pos = sites_list(modified(i),2:3);
    
    [nearest,~] = neighbor_sites(pos(1),pos(2),sublattice);
    [n_near, ~] = size(nearest);
    
    % Loop through nearest neighbor lattice sites
    for k = 1 : n_near
        
        posr = nearest(k,1);
        posc = nearest(k,2);
        neigh_sublattice = nearest(k,3);
        
        [lia,locb] =   ismember([neigh_sublattice,posr,posc],sites_list(:,[1,2,3]),'rows');
        
        if locb > 0
            if (sites_list(locb,4) == 0) %If the site is unoccupied, rate  will be modified if M is neighbor
                b1_b2 = count_M_neighbors(posr,posc,neigh_sublattice,C);
                if b1_b2 > 0
                    sumM = [M2(3)+M2(2)+M2(1) M2(3)+M2(2) M2(3);
                        M1(3)+M1(2)+M1(1) M1(3)+M1(2) M1(3)];
                    
                    RateTotal = km*F_f^(b1_b2-1)*P^b1_b2*sumM(neigh_sublattice,b1_b2) + ki*P*I(neigh_sublattice);
                    
                    sites_list(locb,5) = RateTotal; %Rate updated
                else
                    sites_list(locb,5) = 0;
                end  
            end
        end
        
        if ((~any(lia)) && C(posr,posc,neigh_sublattice) == 0) % If lattice site is not in sites_list, add it
            
            if i <= 2 + sitesLength %If this site has a M neighbor
                b1_b2 = count_M_neighbors(posr,posc,neigh_sublattice,C);
                
                sumM = [M2(3)+M2(2)+M2(1) M2(3)+M2(2) M2(3);
                    M1(3)+M1(2)+M1(1) M1(3)+M1(2) M1(3)];
                
                RateTotal = km*F_f^(b1_b2 - 1)*P^b1_b2*sumM(neigh_sublattice,b1_b2) + ki*P*I(neigh_sublattice);
                sites_list = [sites_list; neigh_sublattice,posr,posc,0,RateTotal,0];
                
                
            end
        end
        
    end
end
end

function [updated_sites_list,modified_sites_index] = sites_change(C,sublattice,modified_sites,sites_list,site_status,inhibitor,n_Species)
% This function updates the sites_list for the simulation and returns an
% updated list of potential sites for reactions

global P F_f km ki Lx Ly Lz NA kir;

Species = n_Species./(Lx*Ly*Lz*1e3*1e-30*NA);
% Species = [M1_3 M1_2 M1_1 M1_0 M2_3 M2_2 M2_1 M2_0 I_1 I_2 I1I2];
M1data=Species(1:4);
M1 = @(x) M1data(4-x);
M2data=Species(5:8);
M2 = @(x) M2data(4-x);
I=Species(9:11);
I = [I(2) I(1)];

updated_sites_list = sites_list;

% Update properties of added or removed lattice site
% I or M update

[Length_modified,~] = size(modified_sites);

for i = 1 : Length_modified
    pos = [modified_sites(i,1) modified_sites(i,2)];
    
    [lia,locb] = ismember([sublattice,pos(1),pos(2),site_status],updated_sites_list(:,[1,2,3,4]),'rows'); % find occupied site in sites_list
    
    if locb > 0
        original_state = updated_sites_list(locb,4);
        updated_sites_list(locb,4) = ~original_state;
        modified_sites_index = locb;
    end
    
    if any(~lia)
        updated_sites_list = [updated_sites_list; sublattice,pos(1),pos(2),1,kir*P,3];
    end
    
end

pos = [modified_sites(1,1) modified_sites(1,2)];

[~,locb] = ismember([sublattice,pos(1),pos(2),~site_status],updated_sites_list(:,[1,2,3,4]),'rows');

[nearest,~] = neighbor_sites(pos(1),pos(2),sublattice);
[n_near, ~] = size(nearest);

if site_status == 0 % M/I addition
    %Update rate and I/M indicator
    
    if inhibitor < 3
        updated_sites_list(locb,6) = 1;
        updated_sites_list(locb,5) = 0;
    else
        updated_sites_list(locb,6) = 3;
        updated_sites_list(locb,5) = kir*P;
    end
    
    %Loop through the neighbors
    for k = 1 : n_near
        posr = nearest(k,1);
        posc = nearest(k,2);
        neigh_sublattice = nearest(k,3);
        
        % Try to find neighoring site in updated_sites_list
        [lia,locb1] = ismember([neigh_sublattice,posr,posc],updated_sites_list(:,[1,2,3]),'rows');
        
        if (any(~lia) && inhibitor <3 && C(posr,posc,neigh_sublattice) == 0) % M neighbor, rate will be enhanced
            b1_b2 = count_M_neighbors(posr,posc,neigh_sublattice,C);
            
            sumM = [M2(3)+M2(2)+M2(1) M2(3)+M2(2) M2(3);
                M1(3)+M1(2)+M1(1) M1(3)+M1(2) M1(3)];
            
            RateTotal = km*F_f^(b1_b2-1)*P^b1_b2*sumM(neigh_sublattice,b1_b2) + ki*P*I(neigh_sublattice);
            updated_sites_list = [updated_sites_list; neigh_sublattice,posr,posc,0,RateTotal,0];
            
        end
        
        if locb1 > 0
            
            if (updated_sites_list(locb1,4)==0 && updated_sites_list(locb1,6)==0) %Site is unoccupied and it's new neighbor is M
                
                b1_b2 = count_M_neighbors(posr,posc,neigh_sublattice,C);
                
                if b1_b2 >0
                    sumM = [M2(3)+M2(2)+M2(1) M2(3)+M2(2) M2(3);
                        M1(3)+M1(2)+M1(1) M1(3)+M1(2) M1(3)];
                    
                    RateTotal = km*F_f^(b1_b2-1)*P^b1_b2*sumM(neigh_sublattice,b1_b2) + ki*P*I(neigh_sublattice);
                    
                    updated_sites_list(locb1,5) = RateTotal;
                else
                    updated_sites_list(locb1,5) = 0;
                end
                
                
            end
        end
    end
    
    %Running a loop through the inhibitor's neighbors
    if Length_modified > 1
        for i = 2 : Length_modified
            pos = [modified_sites(i,1) modified_sites(i,2)];
            [nearest_n,~] = neighbor_sites(pos(1),pos(2),sublattice);
            [n_nearn,~] = size(nearest_n);
            
            for j = 1:n_nearn
                [~,locb2] =   ismember([nearest_n(j,1),nearest_n(j,2),nearest_n(j,3)],updated_sites_list(:,[1,2,3]),'rows');
                if locb2 > 0
                    
                    b1_b2 = count_M_neighbors(pos(1),pos(2),sublattice,C);
                    if b1_b2 >0
                        sumM = [M2(3)+M2(2)+M2(1) M2(3)+M2(2) M2(3);
                            M1(3)+M1(2)+M1(1) M1(3)+M1(2) M1(3)];
                        
                        RateTotal = km*F_f^(b1_b2-1)*P^b1_b2*sumM(neigh_sublattice,b1_b2) + ki*P*I(neigh_sublattice);
                        
                        updated_sites_list(locb2,5) = RateTotal;
                    else
                        updated_sites_list(locb2,5) = 0;
                    end
                end
            end
        end
    end
    
    
    
else %I removal
    updated_sites_list(locb,6) = 0; %As site is removed
    
    posr = pos(1);
    posc = pos(2);
    neigh_sublattice = updated_sites_list(locb,1);
    
    b1_b2 = count_M_neighbors(posr,posc,neigh_sublattice,C);
    
    if b1_b2 > 0
        sumM = [M2(3)+M2(2)+M2(1) M2(3)+M2(2) M2(3);
            M1(3)+M1(2)+M1(1) M1(3)+M1(2) M1(3)];
        
        RateTotal = km*F_f^(b1_b2-1)*P^b1_b2*sumM(neigh_sublattice,b1_b2) + ki*P*I(neigh_sublattice);
        
        updated_sites_list(locb,5) = RateTotal;
    else
        updated_sites_list(locb,5) = 0;
    end
    
    for k = 1:n_near
        
        posr = nearest(k,1);
        posc = nearest(k,2);
        neigh_sublattice = nearest(k,3);
        
        % Try to find neighoring site in updated_sites_list
        [~,locb3] =   ismember([neigh_sublattice,posr,posc],updated_sites_list(:,[1,2,3]),'rows');
        
        if (locb3 > 0) % if found, update the number of occupied neighbors to the site
            
            b1_b2 = count_M_neighbors(updated_sites_list(locb3,2),updated_sites_list(locb3,3),updated_sites_list(locb3,1),C);
            if updated_sites_list(locb3,5) > 0 && b1_b2 > 0
                
                sumM = [M2(3)+M2(2)+M2(1) M2(3)+M2(2) M2(3);
                    M1(3)+M1(2)+M1(1) M1(3)+M1(2) M1(3)];
                
                RateTotal = km*F_f^(b1_b2-1)*P^b1_b2*sumM(neigh_sublattice,b1_b2) + ki*P*I(neigh_sublattice);
                
                updated_sites_list(locb3,5) = RateTotal;
            end
        end
    end
end

end

function y=pbc(x)
% This function implements periodic boundary conditions for lattice
% indices

% Input: row/column position before wrapping
% Output: row/column position after periodic boundary wrapping

global N;
if (x==0 || x==N)
    y=N;
else
    y=mod(x,N);
end
end

function [nearest,next_nearest] = neighbor_sites(i,j,sublattice)
% This function returns the nearest and next-nearest lattice sites for a
% given site (i, j, sublattice)

% Inputs: row (i), column (j), and sublattice (1 or 2)
% Outputs: list of nearest neighbor and next-nearest neighbor lattice sites
% in the format (row, column, sublattice)

if (sublattice==1)
    
    nearest = [pbc(i-1),    pbc(j-1+(~mod(i,2))),        2;
        pbc(i-1),    pbc(j+(~mod(i,2))),          2;
        i,           j,                           2];
    
    next_nearest = [pbc(i-1),    pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i-1),                pbc(j+(~mod(i,2))),          sublattice;
        i,                       pbc(j-1),                    sublattice;
        i,                       pbc(j+1),                    sublattice;
        pbc(i+1),                pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i+1),                pbc(j+(~mod(i,2))),          sublattice];
    
elseif (sublattice==2)
    
    nearest = [i,       j,                       1;
        pbc(i+1),       pbc(j+(~mod(i,2))),      1;
        pbc(i+1),       pbc(j-1+(~mod(i,2))),    1];
    
    next_nearest = [pbc(i-1),       pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i-1),                   pbc(j+(~mod(i,2))),          sublattice;
        i,                          pbc(j-1),                    sublattice;
        i,                          pbc(j+1),                    sublattice;
        pbc(i+1),                   pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i+1),                   pbc(j+(~mod(i,2))),          sublattice];
end
end

function y=count_M_neighbors(i,j,species,C)

%This function gives the number of M occupied neighbors for a particular
%site

[nearest,~] = neighbor_sites(i,j,species);
[n_near, ~] = size(nearest);

y=0;
for k=1:n_near
    posr = nearest(k,1);
    posc = nearest(k,2);
    neigh_spec = nearest(k,3);
    if (C(posr,posc,neigh_spec)==neigh_spec)
        y=y+1;
    end
end

end

function y = count_neighbors(i,j,species,C)

%This function gives the number of M occupied neighbors for a particular
%site

[nearest,~] = neighbor_sites(i,j,species);
[n_near, ~] = size(nearest);

y=0;
for k=1:n_near
    posr = nearest(k,1);
    posc = nearest(k,2);
    neigh_spec = nearest(k,3);
    if C(posr,posc,neigh_spec)>0
        y=y+1;
    end
end

end

function [occupancy,C0,num_atoms] = init(C0,P1P2_type,num_atoms,pair_sel)

% This function adds the P molecules on the sublattices and inhibitors if
% needed (according to the P1P2_type)

center1_1 = pair_sel(1); center1_2 = pair_sel(2);
center2_1 = pair_sel(3); center2_2 = pair_sel(4);

[nearest1,~] = neighbor_sites(center1_1,center1_2,1);
[n_near1, ~] = size(nearest1);

[nearest2,~] = neighbor_sites(center2_1,center2_2,2);
[n_near2, ~] = size(nearest2);

occupancy = [1 center1_1 center1_2];
C0(center1_1,center1_2,1) = 1;
occupancy = [occupancy; 2 center2_1 center2_2];
C0(center2_1,center2_2,2) = 2;
num_atoms = num_atoms + [1,1];

% Adding the inhibitor neighbors according to the P1P2_type

if (P1P2_type(1) == 0) %M1(1) addition
    for i = 1:n_near1
        if (C0(nearest1(i,1),nearest1(i,2),nearest1(i,3)) == 0)
            C0(nearest1(i,1),nearest1(i,2),nearest1(i,3)) = nearest1(i,3) + 2;
            num_atoms = num_atoms + [0,1];
            occupancy = [occupancy; 2 nearest1(i,1) nearest1(i,2)];
        end
    end
end

if (P1P2_type(2) == 0) %M2(1) addition
    for i = 1:n_near2
        if (C0(nearest2(i,1),nearest2(i,2),nearest2(i,3)) == 0)
            C0(nearest2(i,1),nearest2(i,2),nearest2(i,3)) = nearest2(i,3) + 2;
            num_atoms = num_atoms + [1,0];
            occupancy = [occupancy; 1 nearest2(i,1) nearest2(i,2)];
        end
    end
end

if (P1P2_type(1) == 1) %M1(2) addition
    for i = 1:n_near1
        if (C0(nearest1(i,1),nearest1(i,2),nearest1(i,3)) == 0)
            C0(nearest1(i,1),nearest1(i,2),nearest1(i,3)) = nearest1(i,3) + 2;
            num_atoms = num_atoms + [0,1];
            occupancy = [occupancy; 2 nearest1(i,1) nearest1(i,2)];
            break;
        end
    end
end

if (P1P2_type(2) == 1) %M2(2) addition
    for i = 1:n_near2
        if (C0(nearest2(i,1),nearest2(i,2),nearest2(i,3)) == 0)
            C0(nearest2(i,1),nearest2(i,2),nearest2(i,3)) = nearest2(i,3) + 2;
            num_atoms = num_atoms + [1,0];
            occupancy = [occupancy; 1 nearest2(i,1) nearest2(i,2)];
            break;
        end
    end
end

% M1(0) and M2(0) is not possible as M1 and M2 are reacting with each other

end

function []=visualize(C,N,z,reac_sel,Species)
% This function saves the current state of the system, given the surface
% state matrix, as an XYZ file
% Inputs: Surface state matrix (C), Size of grid (N), iteration number(z),
% reactions selected matrix (reac_sel), concentration matrix (Species)
% Outputs: .xyz file

global lattice_const basedirectory;

% Define box sizes in the x and y directions

fid_xyz = fopen(      [basedirectory,'poreEnh_18C_final','.xyz'],'a');
fid_txt = fopen(      [basedirectory,'/poreEnh_18C_final','.txt'],'a');

% Create variables for number of total atoms and number of occupied atoms

num_total_atoms = N*N*2;
num_occupied_atoms = sum(sum(C(:,:,1)==1))+sum(sum(C(:,:,1)==3))+sum(sum(C(:,:,2)==2))+sum(sum(C(:,:,2)==4));


fprintf(fid_xyz,[num2str(num_occupied_atoms),'\n Iteration number: ', num2str(z), ' Reaction Selected: '  ...
    ,num2str(reac_sel),  ' \n']);
fprintf(fid_txt,'%d %d %d %d %d %d %d %d %d %d %d %d %d \n',...
    z, reac_sel,Species(1),Species(2),Species(3),Species(4),Species(5),Species(6),...
    Species(7),Species(8),Species(9),Species(10),Species(11));

count_occupied_atoms=0;

% Cycle through all locations on the 2D lattice
for i=1:N
    for j=1:N
      
        if C(i,j,1)==1
            coord_sublattice_1 = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2  ];
            fprintf(fid_xyz,['M ', num2str(coord_sublattice_1(1)), '  ', num2str(coord_sublattice_1(2)), '   0.000\n']);
            
        elseif C(i,j,1)==3
            coord_sublattice_1 = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2  ];
            fprintf(fid_xyz,['I ', num2str(coord_sublattice_1(1)), '  ', num2str(coord_sublattice_1(2)), '   0.000\n']);
            
        end
        

        if C(i,j,2)==2
            coord_sublattice_2 = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];
            fprintf(fid_xyz,['N ', num2str(coord_sublattice_2(1)), '  ', num2str(coord_sublattice_2(2)), '   0.000\n']);
            
        elseif C(i,j,2)==4
            coord_sublattice_2 = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];
            fprintf(fid_xyz,['J ', num2str(coord_sublattice_2(1)), '  ', num2str(coord_sublattice_2(2)), '   0.000\n']);
            

        end
    end
end


fclose('all');
end