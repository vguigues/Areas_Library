

%
%  Author: Vincent Guigues.
%
%  File : angle.m
%
%  Purpose : assocites an angle to a point (x,y) on a circle of center
%            P and radius d 

function [theta]=angle(x,y,P,d)

if (y>=P(2))   
    theta=acos((x-P(1))/d);
else
    theta=2*(pi)-acos((x-P(1))/d);
end
