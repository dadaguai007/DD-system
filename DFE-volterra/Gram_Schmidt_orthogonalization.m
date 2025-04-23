%  Gram-Schmidt orthogonalization
function Gram_Schmidt_orthogonalization
q(1,:)=p(1,:);

for i=2:length(X)

for j=1:i-1
    apha(i,j)=(p(i,:).'*q(:,j))/(q(:,j).'*q(:,j));
    k=apha(i,j)*q(i,j);
    L=L+k;
end

q(i,:)=p(i,:)-L;
Q(i)=(y.'*q(i,:)).^2/N*(q(i,:).'*q(i,:));
end
