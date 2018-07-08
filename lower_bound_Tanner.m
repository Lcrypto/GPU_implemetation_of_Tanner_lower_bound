%To estimate lower bound of code distance and heuristics optimization
%

function [H distanceVN distanceCN second_over_firt_eigen]= lower_bound_Tanner(qcHFileName)

fid = fopen(qcHFileName, 'r');


n = fscanf(fid, '%d', [1 1]); %VN number
m = fscanf(fid, '%d', [1 1]); %CN number
z = fscanf(fid, '%d', [1 1]); 

I = sparse(eye(z));
Z = sparse(zeros(z));
H = sparse([]);

for i = 1:m
    lH = sparse([]);
    for j = 1:n
        shift = fscanf(fid, '%d', [1 1]);
        if shift == -1
            lH = [lH Z];
        else
            lH = [lH circshift(I, [0 shift])];
        end
    end
    H = [H; lH];
    

end

H_=[zeros(z*m,z*m),H; H',zeros(z*n,z*n)];
[V] = eig(H_);
first_eigen=max(V);
second_eigen=max(V(V~=max(V)));
distanceVN=n*(2*n-second_eigen)/(n*m-second_eigen);
distanceCN=2*n*(2*n+n-2-second_eigen)/(n*(n*m-second_eigen));
second_over_firt_eigen=second_eigen/first_eigen;
end


