clear all;
close all;
clc;

% Initial guess for the reaction constants
 k_init=[70	1.2 301 1.2 2535 194];

 options = optimoptions(@lsqnonlin,'Display','iter');
 k_opt = lsqnonlin(@Residual,k_init,zeros(1,6),inf*ones(1,6),options);

t_span = 0:0.001:24;

% Initial Concentrations of the species 
% C_init = [M1_3, M1_2, M1_1, M1_0, M2_3, M2_2, M2_1, M2_0,
% I1, I2, I1I2, P1_2, P1_1, P1_0, P2_2, P2_1, P2_0]

C_init0=[0.005,0,0,0,0.005,0,0,0,0.06,0.06,0,0,0,0,0,0,0];                  % 12C Case
C_init1=[0.005,0,0,0,0.005,0,0,0,0.03,0.03,0,0,0,0,0,0,0];                  % 6C Case
C_init2=[0.005,0,0,0,0.005,0,0,0,0.015,0.015,0,0,0,0,0,0,0];                % 3C Case
C_init3=[0.005,0,0,0,0.005,0,0,0,0,0,0,0,0,0,0,0,0];                        % 0C Case
C_init4=[0.005,0,0,0,0.005,0,0,0,0.12,0.12,0,0,0,0,0,0,0];                  % 24C Case
C_init5=[0.005,0,0,0,0.005,0,0,0,0.075,0.075,0,0,0,0,0,0,0];                % 15C Case
C_init6=[0.005,0,0,0,0.005,0,0,0,0.09,0.09,0,0,0,0,0,0,0];                  % 18C Case
C_init7=[0.005,0,0,0,0.005,0,0,0,0.105,0.105,0,0,0,0,0,0,0];                % 21C Case
C_init8=[0.005,0,0,0,0.005,0,0,0,0.045,0.045,0,0,0,0,0,0,0];                % 9C Case

% Experimental Yield Values

Y_act0 = [0.0 0.029 0.147 0.275 0.377 0.368 0.392 0.471 0.525 0.539 ...
    0.544 0.642 0.711 0.745 0.794 0.853]';                                  % 12C Case
Y_act1 = [0.0 0.447 0.587 0.686 0.738 0.832 0.864 0.893 0.896 0.889]';% 6C Case
Y_act2 = [0.0 0.391 0.552 0.741 0.824 0.860 0.889 0.893 0.943 0.925]';      % 3C Case
Y_act3 = [0.0 0.035 0.060 0.141 0.310 0.409 0.574 0.739 0.774 0.882 0.902...
    0.953 0.966 0.944]';                                                    % 0C Case
Y_act4 = [0.0 0.010 0.051 0.087 0.126 0.131 0.152]';                         % 24C Case
Y_act5 = [0.0 0.070 0.174 0.241 0.328 0.395 0.462 0.503 0.545]';            % 15C Case
Y_act6 = [0.0 0.048 0.137 0.205 0.247 0.324 0.375 0.407]';                  % 18C Case
Y_act7 = [0.0 0.048 0.115 0.131 0.171 0.221 0.270]';                        % 21C Case
Y_act8 = [0.0 0.412 0.550 0.642 0.712 0.770 0.779 0.853 0.848 0.867]';      % 9C Case

% Calculating the concentrations with optimised reaction constants (k_opt) values 
 [t,C0]= ode23(@(t,C) Rate(t,C,k_opt),t_span,C_init0);                  % 12C Case
 [t,C1]= ode23(@(t,C) Rate(t,C,k_opt),t_span,C_init1);                  % 6C Case
 [t,C2]= ode23(@(t,C) Rate(t,C,k_opt),t_span,C_init2);                  % 3C Case
 [t,C3]= ode23(@(t,C) Rate(t,C,k_opt),t_span,C_init3);                  % 0C Case
 [t,C4]= ode23(@(t,C) Rate(t,C,k_opt),t_span,C_init4);                  % 24C Case
[t,C5]= ode23(@(t,C) Rate(t,C,k_opt),t_span,C_init5);                  % 15C Case
[t,C6]= ode23(@(t,C) Rate(t,C,k_opt),t_span,C_init6);                  % 18C Case
[t,C7]= ode23(@(t,C) Rate(t,C,k_opt),t_span,C_init7);                  % 21C Case
[t,C8]= ode23(@(t,C) Rate(t,C,k_opt),t_span,C_init8);                  % 9C Case

% [t,C]=ode23(@kmcsol, t_span, C_init);

% Calculating the yield from the obtained concentrations
% Y_pred = (M1_3(0) - M1_3(t) - M1_2(t) - M1_1(t) - M1_0 (t))/M1_3(0)

 Y_pred0 = (C_init0(1) - C0(:,1) - C0(:,2) - C0(:,3) - C0(:,4))/C_init0(1);          % 12C Case
 Y_pred1 = (C_init1(1) - C1(:,1) - C1(:,2) - C1(:,3) - C1(:,4))/C_init1(1);          % 6C Case
 Y_pred2 = (C_init2(1) - C2(:,1) - C2(:,2) - C2(:,3) - C2(:,4))/C_init2(1);          % 3C Case
 Y_pred3 = (C_init3(1) - C3(:,1) - C3(:,2) - C3(:,3) - C3 (:,4))/C_init3(1);         % 0C Case
 Y_pred4 = (C_init4(1) - C4(:,1) - C4(:,2) - C4(:,3) - C4(:,4))/C_init4(1);          % 24C Case
Y_pred5 = (C_init5(1) - C5(:,1) - C5(:,2) - C5(:,3) - C5(:,4))/C_init5(1);     
Y_pred6 = (C_init6(1) - C6(:,1) - C6(:,2) - C6(:,3) - C6(:,4))/C_init6(1);   
Y_pred7 = (C_init7(1) - C7(:,1) - C7(:,2) - C7(:,3) - C7(:,4))/C_init7(1);  
Y_pred8 = (C_init8(1) - C8(:,1) - C8(:,2) - C8(:,3) - C8(:,4))/C_init8(1);

%Plotting 
 plot(t_span,Y_pred3);                                         % 0C Case
 hold on;
 plot(t_span,Y_pred2);                                        % 3C Case
 plot(t_span,Y_pred1);                                       % 6C Case
 plot(t_span,Y_pred8);
 plot(t_span,Y_pred0);                                       % 12C Case
 plot(t_span,Y_pred5);                                                      %15C Case
 plot(t_span, Y_pred6);                                                     %18C Case
 plot(t_span,Y_pred7);
 plot(t_span, Y_pred4);                                                     %24C Case
 
t_span0 = [0 0.5 1 1.5 2 2.5 3 4 5 6 7 9 11 13 15 24];                      % 12C Case
 t_span1 = [0 0.5 1 2 4 7 10 15 20 24];                                   % 6C Case
 t_span2 = [0 0.5 1 2 4 7 9.5 15 20 24];                                     % 3C Case
 t_span3 = [0 0.0167 0.0332 0.0833 0.167 0.333 0.5 2 4 7 10 15 20 24];                       % 0C Case
t_span4 = [0 4 7 10 15 20 24];                                              % 24C Case
t_span5 = [0 1 2 4 7 10 15 20 24];                                          % 15C Case
t_span6 = [0 2 4 7 10 15 20 24];                                            % 18C Case
t_span7 = [0 4 7 10 15 20 24];                                              % 21C Case
t_span8 = [0 0.5 1 2 4 7 10 15 20 24];                                      % 9C Case
% 

% 
 plot(t_span3,Y_act3,'gs','MarkerSize',10);  % 0C 
 hold on;
 plot(t_span2,Y_act2,'gs','MarkerSize',10);  % 3C
 plot(t_span1,Y_act1,'gs','MarkerSize',10);    % 6C
 plot(t_span8,Y_act8);    % 9C
  plot(t_span0,Y_act0,'gs','MarkerSize',10);      % 12C
  plot(t_span5,Y_act5,'gs','MarkerSize',10);                                                    %15C Case
  plot(t_span6,Y_act6,'gs','MarkerSize',10); %18C Case
  plot(t_span7,Y_act7,'gs','MarkerSize',10); %21C Case
 plot(t_span4,Y_act4,'gs','MarkerSize',10);                                                     % 24C Case
% plot(t_span/60,Y_pred,'gs','MarkerSize',10,'MarkerEdgeColor','blue','MarkerFaceColor','blue');
 xlabel('Time (in hours)');
 ylabel('Yield');
 ylim([0 1]);
 xlim([0 24]);