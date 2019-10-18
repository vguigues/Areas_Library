
%
%  Author: Vincent Guigues.
%
%  File : is_in_same_half_space
%
%  Purpose : Returns 1 if C and P are in the same half space delimited by line (AB)
%            and 0 otherwise

function [b]=is_in_same_half_space(C,P,A,B)

b=0;

if (A(1)==B(1))
    if ((A(1)-P(1))*(A(1)-C(1))>=0)
       b=1;
    end
else
    fP=A(2)+((B(2)-A(2))/(B(1)-A(1)))*(P(1)-A(1));
    fC=A(2)+((B(2)-A(2))/(B(1)-A(1)))*(C(1)-A(1));
    if ((P(2)-fP)*(C(2)-fC)>=0)
       b=1;
    end
end
    