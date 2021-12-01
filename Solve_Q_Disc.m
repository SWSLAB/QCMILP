function [res]=Solve_Q_Disc(n)
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

%% Q Discretization
Q=sdpvar(Ntot,1);
C=sdpvar(Ntot,1);

% Flag Linear Terms
id_sinks=sum(Bd,1)'==1;
id_sources=sum(Bs,1)'==1;
% Q(id_sinks,:)=Qd;     % In this formulation demands flow is a decision, so we need to keep them in the model
C(id_sources,:)=Cs;
Q(Qmin==Qmax)=Qmin(Qmin==Qmax);
C(Cmin==Cmax)=Cmin(Cmin==Cmax);
M=Q.*C;
for i=1:Ntot
    flag_linear(i)=islinear(M(i));
end
num_bilinear=sum(flag_linear==0);
flag_bilinear=~flag_linear;

%% Define the discretization matrix
% In the flow discretization method the grid we start from zero so
% we have intger multiplications 0*delta 1*delta 2*delta... this is importnat to meet water balance
% Also the grid should be unified for all links, otherwise it will not be possible to meet water balance.
% We also inculde the values of the sinks.
X=repmat([linspace(0,max(Qmax),n-length(Qd)) Qd'],num_bilinear,1);

%% Formulation of the MILP
Y=sdpvar(num_bilinear,n,'full');
b=binvar(num_bilinear,n,'full');

Q(flag_bilinear)=(b.*X)*ones(n,1);
M(flag_linear)=Q(flag_linear).*C(flag_linear);
M(flag_bilinear)=Y*ones(n,1);

Obj=f'*Q;
C1=[A*Q==0];
C2=[A*M==0];
C3=[B*C==0];
%C4=[Bs*C==Cs];  % Defined in line C(id_sources,:)=Cs;
C4=[];
C5=[Bd*Q>=Qd];   % In the flow discretization method the quality constraint is relaxed to inequality
C6=[Qmin<=Q<=Qmax];
C7=[Cmin<=C<=Cmax];
C8=[Y>=0];
C9=[Y<=diag(C(flag_bilinear))*X];
C10=[Y<=(diag(Cmax(flag_bilinear))*X).*b];
C11=[Y>=diag(C(flag_bilinear))*X-(diag(Cmax(flag_bilinear))*X).*(1-b)];
C12=[b*ones(n,1)==ones(num_bilinear,1)];
Const=[C1;C2;C3;C4;C5;C6;C7;C8;C9;C10;C11;C12;];

%% Solve MILP
options=sdpsettings('solver','cplex');
options.verbose=1;
sol=solvesdp(Const,Obj,options);

%% Simulate: check the validity of the solution
ns=length(Cs);
nd=length(Qd);

Qsim=double(Q);
Qsim(abs(Qsim)<=1e-6)=0;
Mat_C=[A*diag(Qsim);B;Bs];
RHS_C=zeros(size(Mat_C,1),1);
RHS_C(end-ns+1:end)=Cs;
Csim=pinv(Mat_C)*RHS_C;

res.n=n;
res.Obj=double(Obj);
res.Q=Qsim;
res.C=Csim;
res.err_salinity_eq=max(abs(Mat_C*Csim-RHS_C))/norm(Csim)*100;
res.err_salinity_ineq=max(max(([Cmin-Csim ; Csim-Cmax]),0))/norm(Csim)*100;

Mat_Q=[A;Bd];
RHS_Q=zeros(size(Mat_Q,1),1);
RHS_Q(end-nd+1:end)=Qd;
res.err_flow_eq=max(abs(Mat_Q*Qsim-RHS_Q))/norm(Qsim)*100;
res.err_flow_ineq=max(max(([Qmin-Qsim ; Qsim-Qmax]),0))/norm(Qsim)*100;

res.max_infeasiblity=max([res.err_salinity_eq res.err_salinity_ineq res.err_flow_eq res.err_flow_ineq]);
res.solvertime=sol.solvertime;

if sol.problem==1 % Infeasible problem
    res.Obj=nan;
    res.max_infeasiblity=nan;
end
end
