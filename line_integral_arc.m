
%
%  Author: Vincent Guigues.
%
%  File : line_integral_arc
%
%  Purpose : output L is one half of the line integral of (xdy-ydx) on arc
%            AB where A,B are two points with angles thetaA, thetaB (given
%            by function angle of the library) on a circle of center P and
%            radius R0, whose computation is explained in
%            "Exact computation of the cumulative distribution function of the distance between a 
%            "point and a random variable uniformly distributed in some sets"
%  available on arXiv. 

function [L]=line_integral_arc(A,B,R0,P,thetaA,thetaB)

L=(R0^2)*(thetaB-thetaA)+P(1)*(B(2)-A(2))-P(2)*(B(1)-A(1));
L=0.5*L;
