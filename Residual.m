function err = Residual(k)

%% Experimental Data
% Measurement time
t_span0 = [0 0.5 1 1.5 2 2.5 3 4 5 6 7 9 11 13 15 24];                      % 12C Case
t_span1 = [0 0.5 1 2 4 7 10 13 15 20 24];                                   % 6C Case
t_span2 = [0 0.5 1 2 4 7 9.5 15 20 24];                                     % 3C Case

% Concentration Initialization
C_init0=[0.005,0,0,0,0.005,0,0,0,0.06,0.06,0,0,0,0,0,0,0];                  % 12C Case
C_init1=[0.005,0,0,0,0.005,0,0,0,0.03,0.03,0,0,0,0,0,0,0];                  % 6C Case
C_init2=[0.005,0,0,0,0.005,0,0,0,0.015,0.015,0,0,0,0,0,0,0];                % 3C Case

% Experimental Yield
Y_act0 = [0.0 0.029 0.147 0.275 0.377 0.368 0.392 0.471 0.525 0.539...
    0.544 0.642 0.711 0.745 0.794 0.853]';                                  % 12C Case
Y_act1 = [0.0 0.401 0.587 0.606 0.747 0.824 0.874 0.869 0.833 0.880 0.889]';% 6C Case
Y_act2 = [0.0 0.391 0.552 0.741 0.824 0.860 0.889 0.893 0.943 0.925]';      % 3C Case
Y_act=[Y_act0; Y_act1; Y_act2];

% Calculating the concentrations for each case at the experimental
% measurement time
[t,C0]= ode23(@(t,C) Rate(t,C,k),t_span0,C_init0);                     % 12C Case
[t,C1]= ode23(@(t,C) Rate(t,C,k),t_span1,C_init1);                     % 6C Case
[t,C2]= ode23(@(t,C) Rate(t,C,k),t_span2,C_init2);                     % 3C Case

% Calculating yield from obtained concentrations
% Y_pred = (M1_3(0) - M1_3(t) - M1_2(t) - M1_1(t) - M1_0(t))/M1_3(0)

Y_pred0=(C_init0(1)-C0(:,1)-C0(:,2)-C0(:,3)-C0(:,4))/C_init0(1);            % 12C Case
Y_pred1=(C_init1(1)-C1(:,1)-C1(:,2)-C1(:,3)-C1(:,4))/C_init1(1);            % 6C Case
Y_pred2=(C_init2(1)-C2(:,1)-C2(:,2)-C2(:,3)-C2(:,4))/C_init2(1);            % 3C Case

% Concatenating the yield matrices
Y_pred=[Y_pred0; Y_pred1; Y_pred2];

% Deviation calculation
err=(Y_pred-Y_act);