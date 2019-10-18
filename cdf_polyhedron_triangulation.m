

%
%  Author: Vincent Guigues.
%
%  File : cdf_polyhedron_triangulation
%
%  Purpose : output f is the value of the cumulative distribution function at d of the distance between 
%            P (a vector of size 2) and a random variable with uniform distribution in polygone S with
%            n vertices. triangles are triangles of a triangulation of S.
%            AreaP is the area of P which can obtained running
%            [Crossing_Number,AreaP,dmin,dmax]=polyhedron(S,P,n).
%  See papers [1] "A library to compute the density of the distance between a point and a random
%              variable uniformly distributed in some sets" and
%             [2] "Computation of the cumulative distribution function of the distance between a 
%             "point and a random variable uniformly distributed in some sets"
%  available on arXiv for details.
%  The algorithm from [1] is used to compute f.


function [f]=cdf_polyhedron_triangulation(S,P,d,n,AreaP,triangles)

AreaInter=0;
for i=1:size(triangles,2)
    [Area] = area_intersection_disk_triangle(P,d,[S(triangles(1,i),1),S(triangles(1,i),2)],[S(triangles(2,i),1),S(triangles(2,i),2)],[S(triangles(3,i),1),S(triangles(3,i),2)]);
    AreaInter=AreaInter+Area;
end
f=AreaInter/AreaP;