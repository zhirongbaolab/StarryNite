function d = distanceToPoint(a,b)
% DISTANCE - computes Euclidean distance matrix
%given a a mx1 point and b m x d vector 

s=size(a);

d=(b(1,:)-a(1)).^2;
for i=2:s(1)
    d=d+(b(i,:)-a(i)).^2;
end
d=sqrt(d);


%d=sqrt((b(1,:)-a(1)).^2+(b(2,:)-a(2)).^2);
%d=d+(b(2,:)-a(2)).^2;
%for i=2:s(1)
%    d=d+(b(i,:)-a(i)).^2;
%end
%d=sqrt(d);

    %{

t1=[1;100];
t2=[5;20];
t2=repmat(t2,1,1);
tic
d1=distance(t1,t2);
toc
tic
%d=sqrt((t2(1,:)-t1(1)).^2+(t2(2,:)-t1(2)).^2);
d2=distanceToPoint(t1,t2);
toc


%}