function D = Eu2_distance(a,b)

[d1,n] = size(a);
[d2,m] = size(b);

if d1~=d2
    error('数据维度不同，无法计算')
end
D = zeros(n,m);
for i = 1:m
    D(:,i) = sum(((a-b(:,i)).^2),1);
end
