
%
%  Author: Vincent Guigues.
%
%  File : arealens
%
%  Purpose : output area is the area of one of the lenses given by chord
%            [I1,I2] and the disk of center C and radius "radius".
%  See papers "A library to compute the density of the distance between a point and a random
%              variable uniformly distributed in some sets" and
%             "Computation of the cumulative distribution function of the distance between a 
%             "point and a random variable uniformly distributed in some sets"
%  available on arXiv for details.



function [area] = arealens(C,radius,I1,I2)
I = (I1+I2)/2;
%if ((I(1)==C(1))&&(I(2)==C(2))): theoretical test replaced by the test
%below due to numerical errors that can occur with this test.
if (norm(I-C)<=10^(-8))
    area=0.5*3.14159*radius^2;
else
    W = I-C;
    V = I2-C;
    dot = W*V';
    costheta = dot/(norm(W) * norm(V));
    area = (radius^2)*acos(costheta)-(radius^2)*sqrt(1-(costheta)^2)*costheta;
end