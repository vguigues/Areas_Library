
%
%  Author: Vincent Guigues.
%
%  File : intersections_ball_segment
%
%  Purpose : N is the number of intersections between segment [A,B] and the circle
%            of center P and radius d with A, B, P vectors of size 2. 
%            If N=0, Intersections=[].
%            If $N=1$ Intersections is a row vector giving the coordinates
%            of the intersection.
%            If N=2, the first row of Intersections gives the coordinates of
%            the first intersection and the second row of Intersections
%            gives the coordinates of the second intersection.
%            
%  The implementation follows the pseudo-code given and explained in 
%  "Exact computation of the cumulative distribution function of the distance between a 
%  "point and a random variable uniformly distributed in some sets"
%  available on arXiv.


function [Intersections,N]=intersections_ball_segment(A,B,P,d)

N=0;
Intersections=[];

Delta=(((A-P)*(B-A)')^2)-((norm(B-A))^2)*((norm(P-A))^2-(d^2));

                                                                                       
if (Delta>=0) 
    if (Delta==0)  
        t1=(-(A-P)*(B-A)')/((norm(B-A))^2);
        I1=A+t1*(B-A);                   
        if ((A-I1)*(B-I1)'<=(10^(-10)))
             N=1;
             Intersections=[I1];
        end  
    else
        t1=(-(A-P)*(B-A)'+sqrt(Delta))/((norm(B-A))^2);
        t2=(-(A-P)*(B-A)'-sqrt(Delta))/((norm(B-A))^2);
        
        I1=A+t1*(B-A);
        I2=A+t2*(B-A);
        
        if ((A-I1)*(B-I1)'<=(10^(-10))) 
            if ((A-I2)*(B-I2)'<=(10^(-10))) 
                N=2;
                Intersections=[I1;I2]; 
            else
                N=1;
                Intersections=[I1];   
            end 
        else
            if ((A-I2)*(B-I2)'<=(10^(-10)))
                N=1;
                Intersections=[I2];
            end
        end  
    end
end


