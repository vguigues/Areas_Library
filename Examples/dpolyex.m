
%
%  Author: Vincent Guigues.
%
%  File : dpolyex.m
%
%  Purpose :   To demonstrate how to compute densities when S is a
%              polyhedron.
%

function [dT1,dT2,dR1,dR2,dP1,dP2]=dpolyex(Nbpoints,algo)

%S is a triangle 
%P outside S

P1=[1,1];
P2=[10,1];
P3=[3,4];
S=[P1;P2;P3;P1];
P=[5,0];
[dT1,timeT1,dmin,dmax]=density_polyhedron(S,P,Nbpoints,algo);
step=(dmax-dmin)/Nbpoints;
subplot(3,4,1);
plot(S(:,1),S(:,2),'r-','Linewidth',2);
hold on
plot(P(1),P(2),'o')
text(P(1)+0.9,P(2),'P');
xlim([0 11]);
ylim([-0.5 5]);
subplot(3,4,2);
plot([dmin+step:step:dmax],dT1);

%S is a triangle 
%P inside S

P1=[1,1];
P2=[10,1];
P3=[3,4];
S=[P1;P2;P3;P1];
P=[4,2];
[dT2,timeT2,dmin,dmax]=density_polyhedron(S,P,Nbpoints,algo);
step=(dmax-dmin)/Nbpoints;
subplot(3,4,3);
plot(S(:,1),S(:,2),'r-','Linewidth',2);
hold on
plot(P(1),P(2),'o')
text(P(1)+0.9,P(2),'P');
xlim([0 11]);
ylim([-0.5 5]);
subplot(3,4,4);
%plot([dmin+step:step:dmax],dT2);
plot([dmin+step:step:dmax-step],dT2(1:length(dT2)-1));

%S is a rectangle 
%P outside S

P1=[3,3];
P2=[12,3];
P3=[12,7];
P4=[3,7];
S=[P1;P2;P3;P4;P1];
P=[1,1];
[dR1,timeR1,dmin,dmax]=density_polyhedron(S,P,Nbpoints,algo);
step=(dmax-dmin)/Nbpoints;
subplot(3,4,5);
plot(S(:,1),S(:,2),'r-','Linewidth',2);
hold on
plot(P(1),P(2),'o')
text(P(1)+0.9,P(2),'P');
xlim([0 13]);
ylim([0 8]);
subplot(3,4,6);
plot([dmin+step:step:dmax],dR1);

%S is a rectangle 
%P inside S

P1=[3,3];
P2=[12,3];
P3=[12,7];
P4=[3,7];
S=[P1;P2;P3;P4;P1];
P=[6,5];
[dR2,timeR2,dmin,dmax]=density_polyhedron(S,P,Nbpoints,algo);
step=(dmax-dmin)/Nbpoints;
subplot(3,4,7);
plot(S(:,1),S(:,2),'r-','Linewidth',2);
hold on
plot(P(1),P(2),'o')
text(P(1)+0.9,P(2),'P');
xlim([0 13]);
ylim([0 8]);
subplot(3,4,8);
%plot([dmin+step:step:dmax],dR2);
plot([dmin+step:step:dmax-step],dT2(1:length(dR2)-1));

%A polyhedron S
%P outside S

P1=[1,1];
P2=[3,1];
P3=[5,2];
P4=[7,1];
P5=[8,3];
P6=[6,3];
P7=[7,6];
P8=[4,5];
P9=[1,3];
P10=[2,2];
S=[P1;P2;P3;P4;P5;P6;P7;P8;P9;P10;P1];
P=[4,0];
[dP1,timeP1,dmin,dmax]=density_polyhedron(S,P,Nbpoints,algo);
step=(dmax-dmin)/Nbpoints;
subplot(3,4,9);
plot(S(:,1),S(:,2),'r-','Linewidth',2);
hold on
plot(P(1),P(2),'o')
text(P(1)+0.7,P(2),'P');
xlim([0 9]);
ylim([-0.5 7]);
subplot(3,4,10);
%plot([dmin+step:step:dmax],dP1);
plot([dmin+step:step:dmax-step],dT2(1:length(dP1)-1));

%A polyhedron S
%P inside S

P1=[1,1];
P2=[3,1];
P3=[5,2];
P4=[7,1];
P5=[8,3];
P6=[6,3];
P7=[7,6];
P8=[4,5];
P9=[1,3];
P10=[2,2];
S=[P1;P2;P3;P4;P5;P6;P7;P8;P9;P10;P1];
P=[4,3];
[dP2,timeP2,dmin,dmax]=density_polyhedron(S,P,Nbpoints,algo);
step=(dmax-dmin)/Nbpoints;
subplot(3,4,11);
plot(S(:,1),S(:,2),'r-','Linewidth',2);
hold on
plot(P(1),P(2),'o')
text(P(1)+0.7,P(2),'P');
xlim([0 9]);
ylim([-0.5 7]);
subplot(3,4,12);
%plot([dmin+step:step:dmax],dP2);
plot([dmin+step:step:dmax-step],dT2(1:length(dP2)-1));

