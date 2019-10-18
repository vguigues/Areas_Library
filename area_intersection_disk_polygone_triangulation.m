
%
%  Author: Vincent Guigues.
%
%  File : area_intersection_disk_polygone_triangulation
%
%  Purpose : output AreaInter is the area of the intersection of the disk of center
%            P and radius d with polygone S. n is the number of vertices of S.
%  See paper  "A library to compute the density of the distance between a point and a random
%              variable uniformly distributed in some sets" 
%  available on arXiv for details.
%  The algorithm used to compute the intersection area is described in Section 2 of the paper above. 

function [AreaInter]=area_intersection_disk_polygone_triangulation(S,P,d,n)

triangles = polygon_triangulate (n,S(:,1),S(:,2));

AreaTr =0;
for i=1:size(triangles,2)
    [Area] = area_intersection_disk_triangle(P,d,[S(triangles(1,i),1),S(triangles(1,i),2)],[S(triangles(2,i),1),S(triangles(2,i),2)],[S(triangles(3,i),1),S(triangles(3,i),2)]);
    AreaTr=AreaTr+Area;
end
         
