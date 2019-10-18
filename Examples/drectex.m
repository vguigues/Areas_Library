

%
%  Author: Vincent Guigues.
%
%  File : drectex.m
%
%  Purpose :   To test function density_polyhedron
%              using Green and Triangulation (its 2 variants)
%              algorithms when S is a rectangle and compare
%              the corresponding densities with the true one.
%              The rectangle has side lenghts L and alpha*L with 0<alpha<=1
%              P is the center of the rectangle

function [dG,dT1,d,ErrG,ErrT1]=drectex(Nbpoints,alpha,L)

P1=[0,0];
P2=[L,0];
P3=[L,alpha*L];
P4=[0,alpha*L];
n=4;
P=[L/2,alpha*L/2];
S=[P1;P2;P3;P4;P1];

subplot(1,4,1);
plot(S(:,1),S(:,2),'r-','Linewidth',2);
hold on
plot(P(1),P(2),'o')
text(P(1)+0.3,P(2),'P');
xlim([0 L+1]);
ylim([0 alpha*L+0.5]);


[dG,timeg,dmin,dmax]=density_polyhedron(S,P,Nbpoints,'g');
step=(dmax-dmin)/Nbpoints;
subplot(1,4,2);
plot([dmin+step:step:dmax],dG);
legend('Green');

[dT1,timet1,dmin,dmax]=density_polyhedron(S,P,Nbpoints,'t1');
step=(dmax-dmin)/Nbpoints;
subplot(1,4,3);
plot([dmin+step:step:dmax],dT1);
legend('Triangulation 1');

%Analytic
vpi=3.14159;
Abs=[dmin+step:step:dmax];
for i=1:Nbpoints
    if (Abs(i)<=0.5*alpha*L)
       d(i)=2*vpi*Abs(i)/(alpha*L*L);
    elseif (Abs(i)<=0.5*L)
       d(i)=(1/(alpha*L*L))*(2*vpi*Abs(i)-4*Abs(i)*acos(alpha*L/(2*Abs(i))));
    elseif (Abs(i)<=0.5*L*sqrt(1+alpha^2))
       d(i)=(1/(alpha*L*L))*(2*vpi*Abs(i)-4*Abs(i)*acos(alpha*L/(2*Abs(i)))-4*Abs(i)*acos(L/(2*Abs(i))));
    else
       d(i)=0;
    end
end
subplot(1,4,4);
plot(Abs,d);
legend('Analytic');

ErrG=max(abs(d-dG));
ErrT1=max(abs(d-dT1));
