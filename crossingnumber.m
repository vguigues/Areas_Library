
%
%  Author: Vincent Guigues.
%
%  File : crossingnumber
%
%  Purpose : output Crossing_Number is the crossing number for vector
%            P of size 2 and polygone S with P in the plane containing S.
%            n is the number of vertices of S.
%            dmin is the minimal distance between P and the border of S.
%            The implementation follows the pseudo-code given and explained in 
%            "Exact computation of the cumulative distribution function of the distance between a 
%            "point and a random variable uniformly distributed in some sets"
%            available on arXiv.

function [Crossing_Number,dmin]=crossingnumber(S,P,n)
    
Crossing_Number=0;
dmin=inf; 

for i=1:n
    if (  ( (S(i,2)<=P(2)) && (S(i+1,2)>P(2))  ) || (  (S(i,2)>P(2))   && (S(i+1,2)<=P(2)) ) ) 
        xI=S(i,1)+((S(i+1,1)-S(i,1))/(S(i+1,2)-S(i,2)))*(P(2)-S(i,2));
        if (xI>=P(1)) 
            Crossing_Number=Crossing_Number+1;            
        end
    end
    
    Vec=S(i+1,:)-S(i,:);
    BarP=S(i,:)+ (((P-S(i,:))*Vec')/(norm(Vec)^2))*Vec;
    if ((S(i,:)-BarP)*(S(i+1,:)-BarP)'<=0)
        dmin=min(dmin, norm(BarP-P));
    else
        dmin=min(dmin, min(norm(S(i,:)-P), norm(S(i+1,:)-P)));
    end
    
end

