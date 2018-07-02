%Copyright(c) 2012, USAtyuk Vasiliy 
%All rights reserved.
%
%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are met :
%*Redistributions of source code must retain the above copyright
%notice, this list of conditions and the following disclaimer.
%* Redistributions in binary form must reproduce the above copyright
%notice, this list of conditions and the following disclaimer in the
%documentation and / or other materials provided with the distribution.
%* Neither the name of the <organization> nor the
%names of its contributors may be used to endorse or promote products
%derived from this software without specific prior written permission.
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED.IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
%DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% Application generate QC-LDPC parity-check matrix based on RS-codewords
%
%Tanner R. M. (2001) ‘Minimum distance bounds by graph analysis’, IEEE-IT, Vol. 47,
%pp. 808-821.
%To estimate lower bound of code distance and heuristics optimization second_eigen/first_eigen  to choice best code
% in ensemble, very fast sieving from code distance, especialy important for low rate codes
% as Vontobel proof it tight bound to Trapping Sets pseudocodewords weight too
function [H distanceVN distanceCN second_over_first_eigen]= lower_bound_Tanner_GPU(qcHFileName)

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
G = gpuArray(full(H_));
[V] = eig(G);
first_eigen=max(V);
second_eigen=max(V(V~=max(V)));
distanceVN=n*(2*n-second_eigen)/(n*m-second_eigen);
distanceCN=2*n*(2*n+n-2-second_eigen)/(n*(n*m-second_eigen));
second_over_first_eigen=second_eigen/first_eigen;
fclose(fid);
end


