function D=OffDiag(M, nonan)


if size(M,1) == size(M,2) 
    issym = max(max(M-M')) < 1e-10;
%     issym = sum(sum(M==M'))==numel(M);
else
    issym = 0;
end

if ~issym
    take=eye(size(M))==0;
else
    take=triu(ones(size(M)),1)==1; 
end

D=M(take);

if nargin==2 && nonan==true
    D=D(~isnan(D));
end
 