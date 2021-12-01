function [C]=build_equal_out_matrix(in)
[m n]=size(in);
C=[];
A=(in==-1);
var=1:n;
for i=1:m
    ix=var(A(i,:));
    dim=length(ix);
    D0=zeros(dim,n);
for j=1:dim
    D0(j,ix(j))=1;
end
if size(D0,1)>1
B=diff(D0);
else
if    isempty(ix) | dim==1
    B=[];
else
B=D0;
end
end
C=[C ;B];
end
