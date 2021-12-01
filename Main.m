%% ####################################################################################################################
% Code for the paper:
% Mixed-Integer Linear Programs for Optimizing Multi-Source Water Supply Systems
% By Mashor Housh, PhD
% University of Haifa, mhoush@univ.haifa.ac.il
%% ####################################################################################################################
% This code requires:
% YALMIP toolbox: https://yalmip.github.io/
% Optitoolbox: https://inverseproblem.co.nz/OPTI/
% Developed under Matlab 2018b
% ####################################################################################################################
clc
clear
Nmax=800;
%% Solve NLP
res_NLP=Solve_NLP();

%% Solve AC-Discretization
figure
vec_input=20:20:Nmax;
for i=1:length(vec_input)
    res_AC(i)=Solve_C_Disc(vec_input(i));
    subplot(2,1,1)
    rel_err=abs([res_AC(1:i).Obj]-res_NLP.Obj)/res_NLP.Obj*100;
    plot([res_AC(1:i).n],rel_err)
    set(gca, 'YScale', 'log')
    xlim([5,500])
    subplot(2,1,2)
    plot([res_AC(1:i).n],[res_AC(1:i).max_infeasiblity])
    set(gca, 'YScale', 'log')
    xlim([5,500])
    drawnow
    pause(1)
end

%% Solve GC-Discretization
figure
vec_input=20:20:Nmax;
Epsilons=(300/20).^(1./(vec_input-14))-1; % Csmax=300, Csmin=20, Nsources=14
for i=1:length(vec_input)
    res_GC(i)=Solve_C_Disc(Epsilons(i));
    subplot(2,1,1)
    rel_err=abs([res_GC(1:i).Obj]-res_NLP.Obj)/res_NLP.Obj*100;
    plot([res_GC(1:i).n],rel_err)
    set(gca, 'YScale', 'log')
    xlim([5,500])
    subplot(2,1,2)
    plot([res_GC(1:i).n],[res_GC(1:i).max_infeasiblity])
    set(gca, 'YScale', 'log')
    xlim([5,500])
    drawnow
    pause(1)
end

%% Solve Q-Discretization
figure
vec_input=20:20:Nmax;
for i=1:length(vec_input)
    res_Q(i)=Solve_Q_Disc(vec_input(i));
    subplot(2,1,1)
    rel_err=abs([res_Q(1:i).Obj]-res_NLP.Obj)/res_NLP.Obj*100;
    plot([res_Q(1:i).n],rel_err)
    set(gca, 'YScale', 'log')
    subplot(2,1,2)
    plot([res_Q(1:i).n],[res_Q(1:i).max_infeasiblity])
    set(gca, 'YScale', 'log')
    drawnow
    pause(1)
end

%% Solve_PLA_SOS2
figure
vec_input=10:10:Nmax/2;
for i=1:length(vec_input)
    res_PLA_SOS2(i)=Solve_PLA_SOS2(vec_input(i),vec_input(i));
    subplot(2,1,1)
    rel_err=abs([res_PLA_SOS2(1:i).Obj]-res_NLP.Obj)/res_NLP.Obj*100;
    plot([res_PLA_SOS2(1:i).Nd1],rel_err)
    set(gca, 'YScale', 'log')
    subplot(2,1,2)
    plot([res_PLA_SOS2(1:i).Nd1],[res_PLA_SOS2(1:i).max_infeasiblity])
    set(gca, 'YScale', 'log')
    drawnow
    pause(1)
end

%% Save Results
clearvars -except res_NLP res_AC res_GC res_Q res_PLA_SOS2
save('Results1.mat')