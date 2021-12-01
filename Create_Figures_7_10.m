%% ####################################################################################################################
% Code for the paper:
% Mixed-Integer Linear Programs for Optimizing Multi-Source Water Supply Systems
% By Mashor Housh, PhD
% University of Haifa, mhoush@univ.haifa.ac.il
%% ####################################################################################################################

clear
clc
close all
load('Results.mat')

%% AC-Discritization
figure
subplot(3,1,1)
rel_err=abs([res_AC.Obj]-res_NLP.Obj)/res_NLP.Obj*100;
plot([res_AC.n],rel_err)
set(gca, 'YScale', 'log')
xlim([20 800])
set(gca,'XTickLabel','','YMinorTick','on')
ylabel('Error (%)');
ylim([0.1 20])
set(gca,'YTick',[0.1 1 10]);

subplot(3,1,2)
plot([res_AC.n],[res_AC.max_infeasiblity])
set(gca, 'YScale', 'log')
xlim([20 800])
set(gca,'XTickLabel','','YMinorTick','on')
ylabel('Max Infeasibility (%)');
ylim([1e-12 1e-10])
set(gca,'YTick',[10^-12 10^-11 10^-10]);
set(gca,'YTickLabel',{'10^{-12}' '10^{-11}' '10^{-10}'},'YMinorTick','on')

subplot(3,1,3)
plot([res_AC.n],[res_AC.solvertime])
set(gca, 'YScale', 'log')
xlim([20 800])
xlabel('$N^C$','Interpreter','latex');
ylabel('Solver Time (sec)');
ylim([0.1 1000])
set(gca,'YTick',[1 10 100 1000]);

%% GC-Discritization
figure
subplot(3,1,1)
rel_err=abs([res_GC.Obj]-res_NLP.Obj)/res_NLP.Obj*100;
plot([res_GC.n],rel_err)
set(gca, 'YScale', 'log')
xlim([20 800])
set(gca,'XTickLabel','','YMinorTick','on')
ylabel('Error (%)');
ylim([0.1 20])
set(gca,'YTick',[0.1 1 10]);

subplot(3,1,2)
plot([res_GC.n],[res_GC.max_infeasiblity])
set(gca, 'YScale', 'log')
xlim([20 800])
set(gca,'XTickLabel','','YMinorTick','on')
ylabel('Max Infeasibility (%)');
ylim([1e-12 1e-10])
set(gca,'YTick',[10^-12 10^-11 10^-10]);
set(gca,'YTickLabel',{'10^{-12}' '10^{-11}' '10^{-10}'},'YMinorTick','on')

subplot(3,1,3)
plot([res_GC.n],[res_GC.solvertime])
set(gca, 'YScale', 'log')
xlim([20 800])
xlabel('$N^C$','Interpreter','latex');
ylabel('Solver Time (sec)');
ylim([0.1 1000])
set(gca,'YTick',[1 10 100 1000]);

%% Q-Discritization
figure
subplot(3,1,1)
rel_err=abs([res_Q.Obj]-res_NLP.Obj)/res_NLP.Obj*100;
plot([res_Q.n],rel_err)
set(gca, 'YScale', 'log')
xlim([20 800])
set(gca,'XTickLabel','','YMinorTick','on')
ylabel('Error (%)');
ylim([0.1 20])
set(gca,'YTick',[0.1 1 10 ]);

subplot(3,1,2)
plot([res_Q.n],[res_Q.max_infeasiblity])
set(gca, 'YScale', 'log')
xlim([20 800])
set(gca,'XTickLabel','','YMinorTick','on')
ylabel('Max Infeasibility (%)');
ylim([0.01 3])
set(gca,'YTick',[0.01 0.1 1]);

subplot(3,1,3)
plot([res_Q.n],[res_Q.solvertime])
set(gca, 'YScale', 'log')
xlim([20 800])
xlabel('$N^Q$','Interpreter','latex');
ylabel('Solver Time (sec)');
ylim([0.1 1000])
set(gca,'YTick',[1 10 100 1000]);

%% PLA-Discritization
figure
subplot(3,1,1)
rel_err=abs([res_PLA_SOS2.Obj]-res_NLP.Obj)/res_NLP.Obj*100;
plot([res_PLA_SOS2.Nd1]*2,rel_err)
set(gca, 'YScale', 'log')
xlim([20 800])
set(gca,'XTickLabel','','YMinorTick','on')
ylabel('Error (%)');
ylim([0.00001 1])
set(gca,'YTick',[0.00001 0.001 0.1 ]);

subplot(3,1,2)
plot([res_PLA_SOS2.Nd1]*2,[res_PLA_SOS2.max_infeasiblity])
set(gca, 'YScale', 'log')
xlim([20 800])
set(gca,'XTickLabel','','YMinorTick','on')
ylabel('Max Infeasibility (%)');

ylim([1e-12 1e-1])
set(gca,'YTick',[10^-12 10^-6 10^-1]);
set(gca,'YTickLabel',{'10^{-12}' '10^{-6}' '10^{-1}'},'YMinorTick','on')

subplot(3,1,3)
plot([res_PLA_SOS2.Nd1]*2,[res_PLA_SOS2.solvertime])
set(gca, 'YScale', 'log')
xlim([20 800])
xlabel('$2\cdot{N^v}$','Interpreter','latex');
ylabel('Solver Time (sec)');
ylim([0.1 1000])
set(gca,'YTick',[1 10 100 1000]);