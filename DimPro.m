a= [0 -1 1 0 0 0 1 0 0 0;...
    0 1 -1 0 0 0 -1 0 0 0;...
    0 1 -1 0 -1 -1 0 0 0 1;...
    0 -1 1 0 1 1 0 0 0 -1;...
    0 -1 1 1 0 0 0 -1 0 -1;...
    0 1 -1 -1 0 0 0 1 0 1;...
    1 -1 1 1 0 0 0 0 -1 0;...
    -1 1 -1 -1 0 0 0 0 1 0];

% the variable matrix
M = [0.645 0.424 0.061 2.25 0.6 0.128 0.495 0.175 0.022 1.25;...
    0.00025 0.00025 0.0007 0.0015 0.003 0.00125 0.0015 0.00175 0.0015 0.00075];
X = M';

% the funtional requirements
F = [-0.125 0.135 -0.875 0.899 -0.450 0.470 -2.500 2.516];
b = F';

% calculate the beta
w33 = a.^2*X(:,2).^2;
denominaor = w33.^(0.5);
beta = (b + a(:,:)*X(:,1))./denominaor

% calculate alpha
as = [];
X_temp = X(:,2);
for k = 1:10
    as_temp = a(:,k).*X_temp(k);
    as = [as as_temp];
end    
alpha = as./denominaor

%calculate covariance "rou"
rou = alpha*alpha'

% calculate the probability of Fi
Pro_F = normcdf(-beta,0,1)

% calculate joint probability
Pro_j = zeros(8);
mu = [0,0];
minusbeta = -1.*beta;
Pro_j(1,1)= Pro_F(1);
for i = 3:8
    Pro_j(i,1) = mvncdf([minusbeta(i),minusbeta(1)],mu,[1 rou(i,1);rou(i,1) 1]);
end
for j = 2:7
    Pro_j(j,j)= Pro_F(j);
    Pro_j(8,8)= Pro_F(8);
    for i = j+1:8
        if rou(i,j)== -1.0000
            Pro_j(i,j)= 0;
        else
            Pro_j(i,j) = mvncdf([minusbeta(i),minusbeta(j)],mu,[1 rou(i,j);rou(i,j) 1]);
        end
    end
end

% for the original order 1,2,3,4,5,6,7,8
Pro_j_change = Pro_j;
for i = 1:8
    Pro_j_change(i,i) = 0; 
end
maxPrF1_F2_j = 0;
for m = 1:8
    maxPrF1_F2_j=  maxPrF1_F2_j + max(Pro_j_change(m,:));
end
PrFi = sum(Pro_F)-maxPrF1_F2_j
PrAA = 1-PrFi

% The heuristic order for sum matrix
C = zeros(2,8);
for i = 1:8
    C(1,i)= i;
    for j = 1:8
    C(2,i) = C(2,i) + Pro_j_change(i,j);
    end
end

% sort with the descend order by the 2th culoum
D = C';
E = sortrows(D,-2);

% calculate and sum up the joint probability of herustic order
maxPrFi_Fj_h_order = 0;
Pro_j_change2 = Pro_j_change+ Pro_j_change';
count = [];
for i = 2:8
    for k = 1:i-1        
        count(k) = Pro_j_change2(E(i,1),E(i-k,1));
    end
    maxPrFi_Fj_h_order = maxPrFi_Fj_h_order + max(count);   
end

% the final probability
PrFi_2 = sum(Pro_F)-maxPrFi_Fj_h_order
PrAA_2 = 1-PrFi_2
