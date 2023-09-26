% This program repeatedly calls the KMC routine and stochastically
% add/remove monomer/inhibitor on lattice sites. 

% Re-initialize the random number generator in MATLAB
rng('shuffle')

% Directory in which the COF growth files will be stored
basedir = '..COF growth/';

    dirname = [basedir,'Case12C'];
    status = mkdir(dirname);
    mkdir([dirname,'/path']);



  % Run the KMC algorithm
  [tf,tknock_timeseries,sites_list,C]=KMC_RPT(basedir);
        
  % Store the COF properties in lists
  tf_list = tf;
  tknock_list = [tknock_list; tknock_timeseries(end)];
    
