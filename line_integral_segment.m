
%
%  Author: Vincent Guigues.
%
%  File : line_integral_segment
%
%  Purpose : output L is one half of the line integral of (xdy-ydx) on
%            segment [A,B] with A,B vectors of size 2.

function [L]=line_integral_segment(A,B)
    
L=0.5*(B(2)*A(1)-B(1)*A(2));
    
