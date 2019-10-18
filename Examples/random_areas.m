

%
%  Author: Vincent Guigues.
%
%  File : random_areas.m
%
%  Purpose :   To test Green and Triangulation algorithms to
%  compute the area of NbSim intersections of disks and polyhedrons
%  comparing their accuracy and CPU time. 

function [Errmax,ErrMoy,TimeGreen,TimeTr1]=random_areas(NbSim)

NbVerticess=[10;25;50;80;100;150;200];
TimeGreen=zeros(NbSim,7);
TimeTr1=zeros(NbSim,7);
ErrMoy=zeros(7,1);
Errmax=-inf*ones(7,1);
Radius=1000;

for j=1:7
    NbVertices=NbVerticess(j);
    for k=1:NbSim
        
        P=[-100+200*rand,-100+200*rand];
        d=50+200*rand;
        [S]=generate_polygone(NbVertices,Radius);
        
        NbV=size(S,1);
        
        tic
        triangles = polygon_triangulate (NbV,S(:,1),S(:,2));
        
        AreaTr=0;
        for i=1:size(triangles,2)
            [Area] = area_intersection_disk_triangle(P,d,[S(triangles(1,i),1),S(triangles(1,i),2)],[S(triangles(2,i),1),S(triangles(2,i),2)],[S(triangles(3,i),1),S(triangles(3,i),2)]);
            AreaTr=AreaTr+Area;
        end
        Aux=toc;
        TimeTr1(k,j)=Aux;
        
        S=[S;[S(1,1),S(1,2)]];
        [Crossing_Number,AreaP,dmin,dmax]=polyhedron(S,P,NbV);
        tic
        AreaG=area_intersection_disk_polygone_green(S,P,d,NbV,Crossing_Number,AreaP);
        Aux=toc;
        TimeGreen(k,j)=Aux;
        
        Errmax(j)=max(abs(AreaG-AreaTr),Errmax(j));
        ErrMoy(j)=ErrMoy(j)+abs(AreaG-AreaTr);
    end
    ErrMoy(j)=ErrMoy(j)/NbSim;
end

