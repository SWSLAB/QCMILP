%% ####################################################################################################################
% Code for the paper:
% Mixed-Integer Linear Programs for Optimizing Multi-Source Water Supply Systems
% By Mashor Housh, PhD
% University of Haifa, mhoush@univ.haifa.ac.il
%% ####################################################################################################################

data=xlsread('Topology.xlsx','data_results','A2:J31');
OD=data(:,4:5);
Qmin=data(:,6);
Qmax=data(:,7);
Cs=data(:,8);
source_flag=~isnan(Cs);
Cs=Cs(source_flag);
Cmax=data(:,9);
Cmax(isnan(Cmax))=max(Cs);
Cmin=0*Cmax+min(Cs);
f=data(:,10);
Nnodes=max(OD(:));
Nlinks=size(OD,1);
Nsources=sum(source_flag);
Ndemands=sum((OD(:,2)==Nnodes));
Qd=Qmin(OD(:,2)==Nnodes);
A=zeros(Nnodes,Nlinks);
for i=1:Nlinks
    A(OD(i,:),i)=[-1 1];
end
A(end,:)=[];
Nnodes=Nnodes-1;
A(OD(source_flag,1),:)=[];
Nnodes=Nnodes-Nsources;
B=build_equal_out_matrix(A);
Bs=zeros(Nsources,Nlinks);
sources_ids=find(source_flag);
for i=1:Nsources
    Bs(i,sources_ids(i))=1;
end
Bd=[zeros(Ndemands,Nlinks-Ndemands) eye(Ndemands)];
Ntot=Nlinks;