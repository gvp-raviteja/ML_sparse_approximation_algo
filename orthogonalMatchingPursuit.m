clc;
clear all;

D = [ 1 0 0 0 0.5 0.6533 0.5 0.2706
       0 1 0 0 0.5 0.2706 -0.5 -0.6533
       0 0 1 0 0.5 -0.2706 -0.5 0.6533
       0 0 0 1 0.5 -0.6533 0.5 -0.2706];
 D1 = [1     0     0     0     0     0     0     0 0.3536    0.4904    0.4619    0.4157    0.3536    0.2778    0.1913    0.0975
     0     1     0     0     0     0     0     0 0.3536    0.4157    0.1913   -0.0975   -0.3536   -0.4904   -0.4619   -0.2778
     0     0     1     0     0     0     0     0 0.3536    0.2778   -0.1913   -0.4904   -0.3536    0.0975    0.4619    0.4157
     0     0     0     1     0     0     0     0  0.3536    0.0975   -0.4619   -0.2778    0.3536    0.4157   -0.1913   -0.4904
     0     0     0     0     1     0     0     0 0.3536   -0.0975   -0.4619    0.2778    0.3536   -0.4157   -0.1913    0.4904
     0     0     0     0     0     1     0     0 0.3536   -0.2778   -0.1913    0.4904   -0.3536   -0.0975    0.4619   -0.4157
     0     0     0     0     0     0     1     0 0.3536   -0.4157    0.1913    0.0975   -0.3536    0.4904   -0.4619    0.2778
     0     0     0     0     0     0     0     1 0.3536   -0.4904    0.4619   -0.4157    0.3536   -0.2778    0.1913   -0.0975];

dim = 8;
k = 2;
indices =  randperm(dim,k);
a=zeros(1,dim);
for i=1:size(indices,2)
    a(indices(i))= randperm(100,1);
end
b= a*D';
%a =[0     0     0     0     0     0     0     0     0    56     0     0     0     0    52     0];
b = b';

%% approximating the coefficient matrix A using OMP algorithm, in the equation A*D = X, given X and D as input
R=b;
An=zeros(1,8);
b_new = zeros(4,1);
i=1;
while sqrt(sum(R.^2)) > 0.001 && i<8
    innerProduct=[];
    for index=1:size(D,2)
        innerProduct(index)=dot(R,D(:,index));
    end
    [c,idx] = max(abs(innerProduct));    
    sum1 = 0;
    if i~=1
        for j=1:i-1
            sum1 = sum1 + (dot(D(:,idx),mu(:,j))/sum(mu(:,j).^2))*mu(:,j);
        end
    end
    mu(:,i) = D(:,idx) - sum1;
    An(idx) = dot(mu(:,i),R)/sum(mu(:,i).^2);
    b_new = b_new + (dot(mu(:,i),R)/sum(mu(:,i).^2))*mu(:,i);
    R = R-(dot(mu(:,i),R)/sum(mu(:,i).^2))*mu(:,i);
    i = i+1;
end

error_a = sqrt(sum((a-An).^2));
error_x = sqrt(sum((b-b_new).^2))
b_new_al = wmpalg('OMP',b,D);
error_x_al = sqrt(sum((b-b_new_al).^2));