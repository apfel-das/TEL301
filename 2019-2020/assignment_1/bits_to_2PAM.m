
function  X = bits_to_2PAM ( b ) 
% X = bits_to_2PAM(b) 
% OUTPUT 
%      Xk: Xk, k=0,...,N-1                                   
% INPUT 
%      b: sequence of bits 
% USAGE:
% Map the input bits as shown below:
%   
%   Inp.           Output.
%   0       ->      +1
%   1       ->      -1
%
% Konstantinos T. Pantelis (the whole ECE Dept. as well has a prototype of this function)


X=zeros(size(b)); %Initialise as zeros with respective length

for k=1:length(b) 
    if(b(k)==0)
        X(k)=1; 
    elseif(b(k)==1)
        X(k)=-1;
    else
        disp('Error') 
        return
    end
end
end
