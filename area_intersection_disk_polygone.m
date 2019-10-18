
%
%  Author: Vincent Guigues.
%
%  File : area_intersection_disk_polygone
%
%  Purpose : output area is the area of the intersection of the disk of center
%            P and radius d with polygone S. n is the number of vertices of S.
%            Crossing_Number and AreaP are obtained calling
%            [Crossing_Number,AreaP,dmin,dmax]=polyhedron(S,P,n).
%  See papers "A library to compute the density of the distance between a point and a random
%              variable uniformly distributed in some sets" and
%             "Exact computation of the cumulative distribution function of the distance between a 
%             "point and a random variable uniformly distributed in some sets"
%  available on arXiv for details.

function [area]=area_intersection_disk_polygone(S,P,d,n,Crossing_Number,AreaP,algo)


if (strcmp(algo,'g')==1)
   [area]=area_intersection_disk_polygone_green(S,P,d,n,Crossing_Number,AreaP);
elseif (strcmp(algo,'t1')==1)
    S=S(1:n,:);
    [AreaInter]=area_intersection_disk_polygone_triangulation(S,P,d,n)
elseif (strcmp(algo,'t2')==1)
    
else
    disp('Unknown value for option algo');
end
