
%
%  Author: Vincent Guigues.
%
%  File : areatriangle
%
%  Purpose : output area is the area of the triangle whose vertices are given by
%            vectors A,B,C of size 2.


function [area] = areatriangle(A, B, C)
    area = 0.5 * abs( A(1) * (B(2)-C(2)) + B(1) * (C(2) - A(2)) + C(1) * (A(2) - B(2)) );
end