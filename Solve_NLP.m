function [res]=Solve_NLP()
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

%% Read data
get_data

%% Nonlinear Formulation
fprintf('Solving NLP... \n')
T=1;
Q=sdpvar(Ntot,T,'full');
C=sdpvar(Ntot,T,'full');

Obj=sum(sum(f.*Q));
C1=[A*Q==0];
C2=[A*(Q.*C)==0];
C3=[B*C==0];
C4=[Bs*C==Cs];
C5=[Bd*Q==Qd];
C6=[Qmin<=Q<=Qmax];
C7=[Cmin<=C<=Cmax];
Const=[C1;C2;C3;C4;C5;C6;C7];

%% Solve using multistart random runs
options=sdpsettings('solver','ipopt');
options.verbose=1;
options.ipopt.max_iter=10000;
options.usex0=1;
rand('seed',1000);
for i=1:50
    assign(Q,unifrnd(Qmin,Qmax));
    assign(C,unifrnd(Cmin,Cmax));
    sol(i)=solvesdp(Const,Obj,options);
    Obj_val(i)=double(Obj);
    Q_val{i}=double(Q);
    C_val{i}=double(C);
    flag(i)=sol(i).problem;
    if flag(i)==0
        fprintf('NLP optimal value is: %f \n', Obj_val(i))
    else
        fprintf('Problem in Solving NLP \n')
    end
end
% Identify best run out of 50
ids=find(flag==0);
[~,idx]=min(Obj_val(ids));
res.Obj=Obj_val(ids(idx));
res.Q=Q_val{ids(idx)};
res.C=C_val{ids(idx)};