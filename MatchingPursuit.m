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

%% generating random input
dim = 8;
k = 3;
indices =  randperm(dim,k);
a=zeros(1,dim);
for i=1:size(indices,2)
    a(indices(i))= randperm(100,1);
end
b= a*D';
%b=[100;50;-50;-100;100;50;-50;-100];
b = b';

%% approximating the coefficient matrix A using MP algotrithm, in the equation A*D = X, given X and D as input
D_org = D;
R=b;
An = zeros(1,8);
b_new = zeros(4,1);
i= 0;
while sqrt(sum(R.^2)) > 1 && i<8
    innerProduct=[];
    for index=1:size(D,2)
        innerProduct(index)=dot(R,D(:,index));
    end
    [c,idx] = max(abs(innerProduct)); 
    An(idx)= An(idx)+innerProduct(idx);
    b_new = b_new + (innerProduct(idx)*D(:,idx));
    R = R-(innerProduct(idx)*D(:,idx));
    %D(:,idx)=0;
    i = i+1;
end

D = D_org;

error_a = sqrt(sum((a-An).^2))
error_x = sqrt(sum((b-b_new).^2))
b_new_al = wmpalg('MP',b,D);
error_x_al = sqrt(sum((b-b_new_al).^2))