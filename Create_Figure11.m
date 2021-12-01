%% ####################################################################################################################
% Code for the paper:
% Mixed-Integer Linear Programs for Optimizing Multi-Source Water Supply Systems
% By Mashor Housh, PhD
% University of Haifa, mhoush@univ.haifa.ac.il
%% ####################################################################################################################

load('Results.mat')
get_data
for i=1:length([res_Q.n])
excess_water(i)=max(res_Q(i).Q(end-3:end)-Qd)
end
figure
plot([res_Q.n],excess_water)
xlim([20 800])
xlabel('$N^Q$','Interpreter','latex');
ylabel('Maximum Excess Water (mcm)');
