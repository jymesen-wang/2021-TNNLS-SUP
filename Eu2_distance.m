function D = Eu2_distance(a,b)

[d1,n] = size(a);
[d2,m] = size(b);

if d1~=d2
    error('����ά�Ȳ�ͬ���޷�����')
end
D = zeros(n,m);
for i = 1:m
    D(:,i) = sum(((a-b(:,i)).^2),1);
end
