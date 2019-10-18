
%
%  Author: Vincent Guigues.
%
%  File : area_intersection_disk_polygone_green
%
%  Purpose : output AreaInter is the area of the intersection of the disk of center
%            P and radius d with polygone S. n is the number of vertices of S.
%            Crossing_Number and AreaP are obtained calling
%            [Crossing_Number,AreaP,dmin,dmax]=polyhedron(S,P,n).
%  See papers [1] "A library to compute the density of the distance between a point and a random
%              variable uniformly distributed in some sets" and
%             [2] "Exact computation of the cumulative distribution function of the distance between a 
%             "point and a random variable uniformly distributed in some sets"
%  available on arXiv for details.
%  The algorithm from [2] is used to compute the intersection area. 



function [AreaInter]=area_intersection_disk_polygone_green(S,P,d,n,Crossing_Number,AreaP)

AreaInter=0;
Intersections=[];
Nb_Intersections=0;
Arcs_Inside=[];

for i=1:n
    if (norm(P-S(i,:))<=d)
        if (norm(P-S(i+1,:))<d)
             [Int]=line_integral_segment(S(i,:),S(i+1,:));
            AreaInter=AreaInter+Int;
        elseif (norm(P-S(i+1,:))>d)
            [Liste,Nb]=intersections_ball_segment(S(i,:),S(i+1,:),P,d);
            Bool=((Liste(1,1)~=S(i,1))||(Liste(1,2)~=S(i,2)));
            if (((Nb==1)&&Bool)||(Nb==2))
                  Pinter=Liste(1,:);
                  [Int]=line_integral_segment(S(i,:),Pinter);
                  Nb_Intersections=Nb_Intersections+1;
                  Intersections=[Intersections;Pinter];
                  Arcs_Inside=[Arcs_Inside;1];
                  AreaInter=AreaInter+Int;
            end
        else
            [Int]=line_integral_segment(S(i,:),S(i+1,:));
            AreaInter=AreaInter+Int;
            if (((S(i+2,1)-S(i+1,1))*(P(1)-S(i+1,1))+(S(i+2,2)-S(i+1,2))*(P(2)-S(i+1,2)))<=0)
                    Nb_Intersections=Nb_Intersections+1;
                    Intersections=[Intersections;S(i+1,:)];
                    Arcs_Inside=[Arcs_Inside;1];
            end
        end
    else
        if (norm(P-S(i+1,:))<d)
            [Liste,Nb]=intersections_ball_segment(S(i,:),S(i+1,:),P,d);
            [Int]=line_integral_segment(Liste(1,:),S(i+1,:));
            Nb_Intersections=Nb_Intersections+1;
            Intersections=[Intersections;Liste(1,:)];
            Arcs_Inside=[Arcs_Inside;0];
            AreaInter=AreaInter+Int;
        elseif (norm(P-S(i+1,:))>d)
            [Liste,Nb]=intersections_ball_segment(S(i,:),S(i+1,:),P,d);
            if (Nb==2)
                [Int]=line_integral_segment(Liste(2,:),Liste(1,:));
                Nb_Intersections=Nb_Intersections+2;
                Intersections=[Intersections;Liste(2,:);Liste(1,:)];
                Arcs_Inside=[Arcs_Inside;0;1];
                AreaInter=AreaInter+Int;
            end
        else
            [Liste,Nb]=intersections_ball_segment(S(i,:),S(i+1,:),P,d);
            if (Nb==1)  
                if (((S(i+2,1)-S(i+1,1))*(P(1)-S(i+1,1))+(S(i+2,2)-S(i+1,2))*(P(2)-S(i+1,2)))>0)
                      Nb_Intersections=Nb_Intersections+1;
                      Intersections=[Intersections;S(i+1,:)];
                      Arcs_Inside=[Arcs_Inside;0];
                end
            elseif (Nb==2)
                [Int]=line_integral_segment(Liste(2,:),Liste(1,:));
                AreaInter=AreaInter+Int;
                Nb_Intersections=Nb_Intersections+1;
                Intersections=[Intersections;Liste(2,:)];
                Arcs_Inside=[Arcs_Inside;0];
                if (((S(i+2,1)-S(i+1,1))*(P(1)-S(i+1,1))+(S(i+2,2)-S(i+1,2))*(P(2)-S(i+1,2)))<=0)
                      Nb_Intersections=Nb_Intersections+1;
                      Intersections=[Intersections;S(i+1,:)];
                      Arcs_Inside=[Arcs_Inside;1];
                end
            end
        end
    end
end

if (Nb_Intersections>0)
    for i=1:Nb_Intersections
        Theta(i)=angle(Intersections(i,1),Intersections(i,2),P,d);
    end
    [SortedTheta,Index]=sort(Theta);
    for i=1:Nb_Intersections
        if (Arcs_Inside(Index(i))==1) 
            A=Intersections(Index(i),:);
            if (i<Nb_Intersections) 
                B=Intersections(Index(i+1),:);
                thetaA=angle(A(1),A(2),P,d);
                thetaB=angle(B(1),B(2),P,d);
                if (thetaA>thetaB) 
                    thetaB=thetaB+2*pi;
                end
                [Int]=line_integral_arc(A,B,d,P,thetaA,thetaB);
            else
                B=Intersections(Index(1),:);
                thetaA=angle(A(1),A(2),P,d);
                thetaB=angle(B(1),B(2),P,d);
                if (thetaA>thetaB) 
                    thetaB=thetaB+2*pi;
                end
                [Int]=line_integral_arc(A,B,d,P,thetaA,thetaB);
            end
            AreaInter=AreaInter+Int;
        end
    end
else
    if (AreaInter~=AreaP) 
        if (mod(Crossing_Number,2)==1) 
            AreaInter=((pi)*(d^2));
        else
            AreaInter=0;
        end
    end
end

