

%
%  Author: Vincent Guigues.
%
%  File : density_polyhedron
%
%  Purpose : output d is a vector of size Np.
%            d(i) is the value of the density of the distance D between 
%            P (a vector of size 2) and a random variable with uniform distribution in polygone S with
%            n vertices. 
%            time is the time needed to compute the density.
%            ]dmin,dmax[ is the support of D.
%            See papers [1] "A library to compute the density of the distance between a point and a random
%                           variable uniformly distributed in some sets" and
%                       [2] "Exact computation of the cumulative distribution function of the distance between a 
%                           "point and a random variable uniformly distributed in some sets"
%            available on arXiv for details.
%            If algo='t1' the algorithm from Section 2 of [1] is used.
%            If algo='g', the algorithm from [2] is used.


function [d,time,dmin,dmax]=density_polyhedron(S,P,Np,algo)

if (strcmp(algo,'g')==1)
    
    d=zeros(1,Np);
    tic
    [n,m]=size(S);
    n=n-1;
    [Crossing_Number,AreaP,dmin,dmax]=polyhedron(S,P,n);
    if (mod(Crossing_Number,2)==1)
        dmin=0;
    end
    FOld=0;
    step=(dmax-dmin)/Np;
    dist=dmin+step;
    
    for  i=1:Np
        [ValCDF]=cdf_polyhedron_green(S,P,dist,n,Crossing_Number,AreaP);
        d(i)=Np*(ValCDF-FOld)/(dmax-dmin);
        dist=dist+step;
        FOld=ValCDF;
    end
    time=toc;
    
elseif (strcmp(algo,'t1')==1)
    
    d=zeros(1,Np);
    tic
    [n,m]=size(S);
    n=n-1;
    [Crossing_Number,AreaP,dmin,dmax]=polyhedron(S,P,n);
    if (mod(Crossing_Number,2)==1)
        dmin=0;
    end
    FOld=0;
    step=(dmax-dmin)/Np;
    dist=[dmin+step:step:dmax];
    S=S(1:n,:);
    triangles = polygon_triangulate (n,S(:,1),S(:,2));
    
    for  i=1:Np
        [ValCDF]=cdf_polyhedron_triangulation(S,P,dist(i),n,AreaP,triangles);
        d(i)=Np*(ValCDF-FOld)/(dmax-dmin);
        FOld=ValCDF;
    end
    time=toc;
    
elseif (strcmp(algo,'t2')==1)
    
else
    disp('Unknown value for option algo');
end
