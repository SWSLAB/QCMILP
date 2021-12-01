%% ####################################################################################################################
% Code for the paper:
% Mixed-Integer Linear Programs for Optimizing Multi-Source Water Supply Systems
% By Mashor Housh, PhD
% University of Haifa, mhoush@univ.haifa.ac.il
%% ####################################################################################################################

figure
q=0:0.1:100;
c=50:0.1:150;
Nv=4;
[qq, cc]=meshgrid(q,c);
subplot(2,2,1)
contourf(qq,cc,qq.*cc,20)
subplot(2,2,2)
x1=25:0.1:125;
f1=x1.^2;
x1a=linspace(25,125,Nv);
f1a=x1a.^2;
plot(x1,f1)
hold all
plot(x1a,f1a,'-o')
subplot(2,2,3)
x2=-25:0.1:75;
f2=-x2.^2;
x2a=linspace(-25,75,Nv);
f2a=-x2a.^2;
plot(x2,f2)
hold all
plot(x2a,f2a,'-o')
subplot(2,2,4)
v1=(cc+qq)/2;
v2=(cc-qq)/2;
for i=1:size(qq,1)
    for j=1:size(qq,2)
        I1(i,j)=interp1(x1a,f1a,v1(i,j));
        I2(i,j)=interp1(x2a,f2a,v2(i,j));
        z(i,j)=I1(i,j)+I2(i,j);
    end
end
err=(z-qq.*cc);
contourf(qq,cc,err,20)