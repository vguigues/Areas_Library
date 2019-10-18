

%
%  Author: Vincent Guigues.
%
%  File : generate_polygone.m
%
%  Purpose :   Generates a nonintersecting polyhedron with at most N vertices
%              points are generated randomly with distances of the vertices 
%              to the origin having uniform distribution on the interval [0,R]
%              [Polygone(i,1);Polygone(i,2)] are the coordinates of ith vertex. 

function [Polygone]=generate_polygone(N,R)

Rs=R*rand(4*N,1);

Thetas=[(pi/2)*rand(N,1);(pi/2)+(pi/2)*rand(N,1);pi+(pi/2)*rand(N,1);1.5*pi+(pi/2)*rand(N,1)];

[Thetas,indexes]=sort(Thetas);
Rs=Rs(indexes);

x=Rs(1)*cos(Thetas(1));
y=Rs(1)*sin(Thetas(1));
Polygone=[x,y];
Current=Thetas(1);
i=2;

while (i<=4*N)
        if (Thetas(i)~=Current)
            x=Rs(i)*cos(Thetas(i));
            y=Rs(i)*sin(Thetas(i));
            Polygone=[Polygone;x,y];
            Current=Thetas(i);
            i=i+1;
        else
            bool=1;
            while ((bool==1)&&(i<=4*N))
                   if (Thetas(i)==Current)
                       i=i+1;
                   else
                       bool=0;
                   end
            end
            if (i<=4*N)        
                x=Rs(i-1)*cos(Thetas(i-1));
                y=Rs(i-1)*sin(Thetas(i-1));
                Polygone=[Polygone;x,y];
                x=Rs(i)*cos(Thetas(i));
                y=Rs(i)*sin(Thetas(i));
                Polygone=[Polygone;x,y];
                Current=Thetas(i);
                i=i+1;
            else
                x=Rs(i-1)*cos(Thetas(i-1));
                y=Rs(i-1)*sin(Thetas(i-1));
                Polygone=[Polygone;x,y]; 
            end
        end
end

Nb=size(Polygone,1);






