

%
%  Author: Vincent Guigues.
%
%  File : polyhedron
%
%  Inputs: S is a polygone in R2 with n vertices and P is a vector of size 2
%          in the plane containing S.
%  Purpose : output Crossing_Number is the crossing number for vector
%            P of size 2 and polygone S.
%            AreaP is the area of P.
%            dmin is the minimal distance between P and the border of S.
%            dmax is the maximal distance between P and the border of S.
%  The implementation follows the pseudo-code given and explained in 
%  "Exact computation of the cumulative distribution function of the distance between a 
%  "point and a random variable uniformly distributed in some sets"
%  available on arXiv.

function [Crossing_Number,AreaP,dmin,dmax]=polyhedron(S,P,n)
    
Crossing_Number=0;
dmax=0;
dmin=inf; 
AreaP=0;

for i=1:n
    [Int]=line_integral_segment(S(i,:),S(i+1,:));
    AreaP=AreaP+Int;
    if (((S(i,2)<P(2)) && (S(i+1,2)>=P(2))) || ((S(i,2)>=P(2))&&(S(i+1,2)<P(2)))) 
        xI=S(i,1)+((S(i+1,1)-S(i,1))/(S(i+1,2)-S(i,2)))*(P(2)-S(i,2));
        if (xI>P(1)) 
            Crossing_Number=Crossing_Number+1;            
        end
    end
    dmax=max(dmax,norm(P-S(i,:)));
    Vec=S(i+1,:)-S(i,:);
    BarP=S(i,:)+ (((P-S(i,:))*Vec')/(norm(Vec)^2))*Vec;
    if ((S(i,:)-BarP)*(S(i+1,:)-BarP)'<=0)
        dmin=min(dmin, norm(BarP-P));
    else
        dmin=min(dmin, min(norm(S(i,:)-P), norm(S(i+1,:)-P)));
    end
end

