function K = F45calc()

% This function calculates the number of free ends of the monomer that is attaching to the
% COF and number of neighbors this monomer will attach to in the COF polymer

% The columns of output K are: 1. Number of neighbors of M1 (x), 2. Number of
% neighbors of P2_1 (b1), 3. Number of neighbors of P2_2 (b2)

K=[]; 

for x = 1 : 3
        
                for b1 = 0 : 3
                    for b2 = 0 : 3
                        
% The first condition checks if the number of COF ends this monomer
% will be connecting to is less than or equal to the free ends of the
% monomer and the second checks if this monomer is connecting to a COF end
% or not

                        if ((b1+b2) <= x && (b1+b2)>0)
                            
                            K = [K ; x b1 b2];
                            
                        end        
                    end
                end
        

end