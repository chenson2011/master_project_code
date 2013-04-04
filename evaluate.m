% Evaluate 

result = [];
for i = 1:2
    
nameA = strcat('network',int2str(i),'.dat');
nameB = strcat('community',int2str(i),'.dat');
A = load(nameA);
B = load(nameB);
cluster = max(B(:,2));



[m,n] = size(A);

matrix = zeros(1000,1000);

for num = 1:m  
    
    matrix(A(num,1),A(num,2)) = A(num,3);
    
end

index = spectral(matrix,cluster);

result = [result,nmi(index',B(:,2)')];

end
result