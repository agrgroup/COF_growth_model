function K = F3calc()

% This function calculates the number of free ends of monomers involved in
% the dimerisation reaction and number of neighbors this dimer will attach
% to in the COF polymer

% Output K has six columns which are: 1. Number of free neighbors of M1 (x), 2. 
% Number of free neighbors of M2 (y), 3. Number of P1_1 neighbors (a1), 4.Number of 
% P1_2 neighbors (a2), 5. Number of P2_1 neighbors (b1), 6. Number of P2_2
% neighbors (b2)

K = []; 

for x = 3 : -1 : 1
    for y = 3: -1 : 1
        for a1 = 0 : y - 1
            for a2 = 0 : y - 1
                for b1 = 0 : x - 1
                    for b2 = 0 : x - 1

% The following condition checks if the number of neighbors in the COF the dimer is
% attaching to is less than the number of free bonds of the M1 and M2

                        if (a1+a2) < y  && (b1+b2) < x 
                            
                            K = [K ; x y a1 a2 b1 b2];
                            
                        end                                 
                    end
                end
            end
        end
    end
end