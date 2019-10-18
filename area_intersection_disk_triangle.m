
%
%  Author: Vincent Guigues.
%
%  File : area_intersection_disk_triangle
%
%  Purpose : output area is the area of the intersection of the disk of center
%            "center" (a row vector of size 2) and radius "radius" with the 
%            triangle whose vertices are given by row vectors A,B,C of size 2. 
%  See paper  "A library to compute the density of the distance between a point and a random
%              variable uniformly distributed in some sets" 
%  available on arXiv for details.
%  The algorithm used to compute the intersection area is described in Section 2 of the paper above. 


function [area] = area_intersection_disk_triangle(center, radius, A, B, C)
    
    [list1, n1] = intersections_ball_segment(A, B, center, radius);
    [list2, n2] = intersections_ball_segment(B, C, center, radius);
    [list3, n3] = intersections_ball_segment(C, A, center, radius);
    code=n1+3*n2+9*n3;
    P = center;
    dA = norm(A-P);
    dB = norm(B-P);
    dC = norm(C-P);
    
    area = 0;
    
    if (code == 0)
        % 000 
        if (dA < radius) && (dB < radius) && (dC < radius)
            area = areatriangle(A,B,C);
        else
            [cn,dmina] = crossingnumber([A;B;C;A], center, 3);
            if ((mod(cn,2)==0)&&(dmina==0))
               cn=1;
            end
            if (mod(cn, 2) == 1)
                area = pi * radius^2;
            end
        end
        
    elseif ((code == 9) || (code == 3) || (code == 1))
        % 001 010 100
         [cn,dmina] = crossingnumber([A;B;C;A], center, 3);
         if ((mod(cn,2)==0)&&(dmina==0))
               cn=1;
         end  
         if mod(cn, 2) == 1
            area = pi * radius^2;
         end
        
    elseif (code == 2)
        % 002
        area = arealens(center, radius, list1(1,:), list1(2,:));
        [cn,dmina] = crossingnumber([A;B;C;A], center, 3);
        if ((mod(cn,2)==0)&&(dmina==0))
             cn=1;
        end   
        if (mod(cn, 2) == 1)
                area = pi * radius^2 - area;
        end
        
    elseif (code == 6)
        % 020
        area = arealens(center, radius, list2(1,:), list2(2,:));
        [cn,dmina] = crossingnumber([A;B;C;A], center, 3);
        if ((mod(cn,2)==0)&&(dmina==0))
             cn=1;
        end  
        if (mod(cn, 2) == 1)
                area = pi * radius^2 - area;
        end
        
    elseif (code == 18)
        % 200
        area = arealens(center, radius, list3(1,:), list3(2,:));
        [cn,dmina] = crossingnumber([A;B;C;A], center, 3);
        if ((mod(cn,2)==0)&&(dmina==0))
             cn=1;
        end  
        if (mod(cn, 2) == 1)
                area = pi * radius^2 - area;
        end
        
    elseif (code == 12)    
        %110
        if ((dA<radius) && (dB<radius) && (dC==radius))
            area = areatriangle(A, B, C);
        elseif ((dA<radius) && (dB<radius) && (dC>radius))
            area= arealens(center, radius, list2(1,:), list3(1,:)) + areatriangle(A, B, list2(1,:))+areatriangle(A, list2(1,:),list3(1,:));
        elseif ((dA>radius) && (dB>radius) && (dC>radius))
            area=pi * radius^2;
        elseif ((dA>radius) && (dB>radius) && (dC<radius))
            [b]=is_in_same_half_space(C,center,list2(1,:), list3(1,:));
            if (b == 0)
               area= pi * radius^2 -arealens(center, radius, list2(1,:), list3(1,:)) + areatriangle(C,list3(1,:) , list2(1,:)); 
            else
               area= arealens(center, radius, list2(1,:), list3(1,:)) + areatriangle(C,list3(1,:) , list2(1,:));
            end
        end
        
    elseif (code == 4)
        %011
        if ((dC<radius) && (dA<radius) && (dB==radius))
            area = areatriangle(A, B, C);
        elseif ((dC<radius) && (dA<radius) && (dB>radius))
            area= arealens(center, radius, list2(1,:), list1(1,:)) + areatriangle(C, A, list1(1,:))+areatriangle(C, list1(1,:),list2(1,:));
        elseif ((dA>radius) && (dC>radius) && (dB>radius))
            area=pi * radius^2;
        elseif ((dA>radius) && (dC>radius) && (dB<radius))
            [b]=is_in_same_half_space(B,center,list2(1,:), list1(1,:));
            if (b == 0)
               area= pi * radius^2 - arealens(center, radius, list2(1,:), list1(1,:)) + areatriangle(B,list1(1,:) , list2(1,:));
            else
               area= arealens(center, radius, list2(1,:), list1(1,:)) + areatriangle(B,list1(1,:) , list2(1,:));
            end
        end
        
    elseif (code == 10)
        %101
        if ((dC<radius) && (dB<radius) && (dA==radius))
            area = areatriangle(A, B, C);
        elseif ((dC<radius) && (dB<radius) && (dA>radius))
            area= arealens(center, radius, list3(1,:), list1(1,:)) + areatriangle(B, C, list3(1,:))+areatriangle(B, list1(1,:),list3(1,:));
        elseif ((dA>radius) && (dC>radius) && (dB>radius))
            area=pi * radius^2;
        elseif ((dB>radius) && (dC>radius) && (dA<radius))
            [b]=is_in_same_half_space(A,center,list1(1,:), list3(1,:));
            if (b == 0)
               area= pi * radius^2 - arealens(center, radius, list1(1,:), list3(1,:)) + areatriangle(A,list1(1,:) , list3(1,:)); 
            else   
               area= arealens(center, radius, list1(1,:), list3(1,:)) + areatriangle(A,list1(1,:) , list3(1,:));
            end
        end
        
    elseif ((code == 11) || (code == 5))
        %201 || 210
        area = arealens(center, radius, list1(1,:), list1(2,:));
        [cn,dmina] = crossingnumber([A;B;C;A], center, 3);
        if ((mod(cn,2)==0)&&(dmina==0))
             cn=1;
        end   
        if (mod(cn, 2) == 1)
                area = pi * radius^2 - area;
        end
    elseif ((code == 21) || (code == 19))
        %012 || 102
        area = arealens(center, radius, list3(1,:), list3(2,:));
        [cn,dmina] = crossingnumber([A;B;C;A], center, 3);
        if ((mod(cn,2)==0)&&(dmina==0))
             cn=1;
        end   
        if (mod(cn, 2) == 1)
                area = pi * radius^2 - area;
        end
    elseif ((code == 7) || (code == 15))
        %120 || 021
        area = arealens(center, radius, list2(1,:), list2(2,:));
        [cn,dmina] = crossingnumber([A;B;C;A], center, 3);
        if ((mod(cn,2)==0)&&(dmina==0))
              cn=1;
        end  
        if (mod(cn, 2) == 1)
                area = pi * radius^2 - area;
        end
    elseif (code == 24)
        %220
        if (norm(C-list3(1,:))<=norm(C-list3(2,:)))
               list31=list3(1,:);
               list32=list3(2,:);
        else
               list31=list3(2,:);
               list32=list3(1,:);
        end
        if (norm(C-list2(1,:))<=norm(C-list2(2,:)))
               list21=list2(1,:);
               list22=list2(2,:);
        else
               list21=list2(2,:);
               list22=list2(1,:);
        end
        
        [b]=is_in_same_half_space(C,center,list32, list22);
        
        if (dC>radius)
           if (b==1)
              area = arealens(center, radius, list31, list21) + arealens(center, radius, list32, list22)+areatriangle(list31, list21, list32) +areatriangle(list21, list32, list22);    
           else
              area = arealens(center, radius, list31, list21) + pi*(radius^2)-arealens(center, radius, list32, list22)+areatriangle(list31, list21, list32) +areatriangle(list21, list32, list22);     
           end           
        else
           if (b==1)
              area = arealens(center, radius, list32, list22)+areatriangle(C, list32, list22);    
           else
              area =  pi*(radius^2)-arealens(center, radius, list32, list22)+areatriangle(C, list32, list22);     
           end 
        end
    elseif (code == 8)
        %022
        if (norm(B-list1(1,:))<=norm(B-list1(2,:)))
               list11=list1(1,:);
               list12=list1(2,:);
        else
               list11=list1(2,:);
               list12=list1(1,:);
        end
        if (norm(B-list2(1,:))<=norm(B-list2(2,:)))
               list21=list2(1,:);
               list22=list2(2,:);
        else
               list21=list2(2,:);
               list22=list2(1,:);
        end
        
        [b]=is_in_same_half_space(B,center,list12, list22);
        
        if (dB>radius)
           if (b==1)
              area = arealens(center, radius, list11, list21) + arealens(center, radius, list12, list22)+areatriangle(list11, list21, list12) +areatriangle(list21, list12, list22);    
           else
              area = arealens(center, radius, list11, list21) + pi*(radius^2)-arealens(center, radius, list12, list22)+areatriangle(list11, list21, list12) +areatriangle(list21, list12, list22);     
           end 
        else
           if (b==1)
              area = arealens(center, radius, list12, list22)+areatriangle(B, list12, list22);    
           else
              area =  pi*(radius^2)-arealens(center, radius, list12, list22)+areatriangle(B, list12, list22);     
           end 
        end
    elseif (code == 20)
        %202
        
        if (norm(A-list1(1,:))<=norm(A-list1(2,:)))
               list11=list1(1,:);
               list12=list1(2,:);
        else
               list11=list1(2,:);
               list12=list1(1,:);
        end
        if (norm(A-list3(1,:))<=norm(A-list3(2,:)))
               list31=list3(1,:);
               list32=list3(2,:);
        else
               list31=list3(2,:);
               list32=list3(1,:);
        end
        
        [b]=is_in_same_half_space(A,center,list12, list32);
        
        if (dA>radius)
           if (b==1)
              area = arealens(center, radius, list11, list31) + arealens(center, radius, list12, list32)+areatriangle(list11, list31, list12) +areatriangle(list31, list12, list32);    
           else
              area = arealens(center, radius, list11, list31) + pi*(radius^2)-arealens(center, radius, list12, list32)+areatriangle(list11, list31, list12) +areatriangle(list31, list12, list32);     
           end 
        else
           if (b==1)
              area = arealens(center, radius, list12, list32)+areatriangle(A, list12, list32);    
           else
              area =  pi*(radius^2)-arealens(center, radius, list12, list32)+areatriangle(A, list12, list32);     
           end     
        end
    elseif (code == 13)
        %111%%%
        if ((dA > radius) && (dB > radius) && (dC > radius)) 
            area = pi*radius^2;
        elseif ((dA==radius) && (dB > radius) && (dC < radius))
            area = arealens(center, radius, list2(1,:), A) + areatriangle(A,C,list2(1,:));
        elseif ((dA==radius) && (dB < radius) && (dC > radius))    
            area = arealens(center, radius, list2(1,:), A) + areatriangle(A,B,list2(1,:));
        elseif ((dB==radius) && (dA > radius) && (dC < radius))
            area = arealens(center, radius, list3(1,:), B) + areatriangle(C,B,list3(1,:));
        elseif ((dB==radius) && (dA < radius) && (dC > radius))    
            area = arealens(center, radius, list3(1,:), B) + areatriangle(A,B,list3(1,:));
        elseif ((dC==radius) && (dB > radius) && (dA < radius))
            area = arealens(center, radius, list1(1,:), C) + areatriangle(A,C,list1(1,:));
        elseif ((dC==radius) && (dB < radius) && (dA > radius))    
            area = arealens(center, radius, list1(1,:), C) + areatriangle(B,C,list1(1,:));
        elseif ((dA<radius) && (dB > radius) && (dC > radius))
            [b]=is_in_same_half_space(A,center,list1(1,:), list3(1,:));
            if (b==1)
                area = arealens(center, radius, list1(1,:), list3(1,:)) + areatriangle(A,list3(1,:),list1(1,:));
            else
                area = pi*(radius^2)-arealens(center, radius, list1(1,:), list3(1,:)) + areatriangle(A,list3(1,:),list1(1,:));
            end
        elseif ((dB<radius) && (dA > radius) && (dC > radius))
            [b]=is_in_same_half_space(B,center,list1(1,:), list2(1,:));
            if (b==1)
                area = arealens(center, radius, list1(1,:), list2(1,:)) + areatriangle(B,list2(1,:),list1(1,:));
            else
                area = pi*(radius^2)-arealens(center, radius, list1(1,:), list2(1,:)) + areatriangle(B,list2(1,:),list1(1,:));
            end
        elseif ((dC<radius) && (dA > radius) && (dB > radius))    
            [b]=is_in_same_half_space(C,center,list2(1,:), list3(1,:));
            if (b==1)
                area = arealens(center, radius, list2(1,:), list3(1,:)) + areatriangle(C,list3(1,:),list2(1,:));
            else
                area = pi*(radius^2)-arealens(center, radius, list2(1,:), list3(1,:)) + areatriangle(C,list3(1,:),list2(1,:));
            end
        end
    elseif (code == 22)
        %211
        if ((dA > radius) && (dB > radius) && (dC > radius))
           [b]=is_in_same_half_space(B,center,A, C); 
           if (b==1)
              area = pi*(radius^2)-arealens(center, radius, list3(1,:), list3(2,:));   
           else
              area = arealens(center, radius, list3(1,:), list3(2,:));    
           end
        elseif ((dA > radius) && (dB < radius) && (dC > radius))
              if (norm(list3(1,:)-A)<=norm(list3(2,:)-A))  
                  l31=list3(1,:);
                  l32=list3(2,:);
              else
                  l31=list3(2,:);
                  l32=list3(1,:);                  
              end
              area = arealens(center, radius, list1(1,:), l31)+arealens(center, radius, list2(1,:), l32)+areatriangle(B,list1(1,:),l31)+areatriangle(B,l32,l31)+areatriangle(B,l32,list2(1,:));
        elseif ((dA==radius) && (dC==radius))
             if (norm(B-center)<radius)
                 area=areatriangle(A,B,C);
             else
                 area=arealens(center, radius, A, C);
             end
        elseif (dC==radius)
             if (norm(B-center)<radius)
                 if (norm(list3(1,:)-C)<=norm(list3(2,:)-C))  
                    l31=list3(2,:);
                 else
                    l31=list3(1,:);                  
                  end
                  area=arealens(center, radius, list1(1,:), l31)+areatriangle(B,list1(1,:), l31)+areatriangle(B,C, l31);
             else
                 [b]=is_in_same_half_space(B,center,A, C);
                 if (b==1)
                     area=pi*(radius^2)-arealens(center,radius,list3(1,:),list3(2,:));
                 else
                     area=arealens(center,radius,list3(1,:),list3(2,:));
                 end
             end
        elseif (dA==radius)
            if (norm(B-center)<radius)
                 if (norm(list3(1,:)-A)<=norm(list3(2,:)-A))  
                    l31=list3(2,:);
                 else
                    l31=list3(1,:);                  
                  end
                  area=arealens(center, radius, list2(1,:), l31)+areatriangle(B,list2(1,:), l31)+areatriangle(B,A, l31);
             else
                 [b]=is_in_same_half_space(B,center,A, C);
                 if (b==1)
                     area=pi*(radius^2)-arealens(center,radius,list3(1,:),list3(2,:));
                 else
                     area=arealens(center,radius,list3(1,:),list3(2,:));
                 end
             end
        end    
    elseif (code == 16)
        %121
        if ((dA > radius) && (dB > radius) && (dC > radius))
           [b]=is_in_same_half_space(A,center,B, C); 
           if (b==1)
              area = pi*(radius^2)-arealens(center, radius, list2(1,:), list2(2,:));   
           else
              area = arealens(center, radius, list2(1,:), list2(2,:));    
           end
        elseif ((dB > radius) && (dA < radius) && (dC > radius))
              if (norm(list2(1,:)-B)<=norm(list2(2,:)-B))  
                  l21=list2(1,:);
                  l22=list2(2,:);
              else
                  l21=list2(2,:);
                  l22=list2(1,:);                  
              end
              area = arealens(center, radius, list1(1,:), l21)+arealens(center, radius, list3(1,:), l22)+areatriangle(A,list1(1,:),l21)+areatriangle(A,l22,l21)+areatriangle(A,l22,list3(1,:)); 
        elseif ((dB==radius) && (dC==radius))
             if (norm(A-center)<radius)
                 area=areatriangle(A,B,C);
             else
                 area=arealens(center, radius, B, C);
             end
        elseif (dC==radius)
             if (norm(A-center)<radius)
                 if (norm(list2(1,:)-C)<=norm(list2(2,:)-C))  
                    l21=list2(2,:);
                 else
                    l21=list2(1,:);                  
                  end
                  area=arealens(center, radius, list1(1,:), l21)+areatriangle(A,list1(1,:), l21)+areatriangle(A,C, l21);
             else
                 [b]=is_in_same_half_space(A,center,B, C);
                 if (b==1)
                     area=pi*(radius^2)-arealens(center,radius,list2(1,:),list2(2,:));
                 else
                     area=arealens(center,radius,list2(1,:),list2(2,:));
                 end
             end
        elseif (dB==radius)
            if (norm(A-center)<radius)
                 if (norm(list2(1,:)-B)<=norm(list2(2,:)-B))  
                    l21=list2(2,:);
                 else
                    l21=list2(1,:);                  
                  end
                  area=arealens(center, radius, list3(1,:), l21)+areatriangle(A,list3(1,:), l21)+areatriangle(A,B, l21);
             else
                 [b]=is_in_same_half_space(A,center,B, C);
                 if (b==1)
                     area=pi*(radius^2)-arealens(center,radius,list2(1,:),list2(2,:));
                 else
                     area=arealens(center,radius,list2(1,:),list2(2,:));
                 end
             end
        end 
    elseif (code == 14)
        %112
        if ((dA > radius) && (dB > radius) && (dC > radius))
           [b]=is_in_same_half_space(C,center,B, A); 
           if (b==1)
              area = pi*(radius^2)-arealens(center, radius, list1(1,:), list1(2,:));   
           else
              area = arealens(center, radius, list1(1,:), list1(2,:));    
           end
        elseif ((dA > radius) && (dC < radius) && (dB > radius))
              if (norm(list1(1,:)-A)<=norm(list1(2,:)-A))  
                  l11=list1(1,:);
                  l12=list1(2,:);
              else
                  l11=list1(2,:);
                  l12=list1(1,:);                  
              end
              area = arealens(center, radius, list3(1,:), l11)+arealens(center, radius, list2(1,:), l12)+areatriangle(C,list3(1,:),l11)+areatriangle(C,l12,l11)+areatriangle(C,l12,list2(1,:)); 
        elseif ((dB==radius) && (dA==radius))
             if (norm(C-center)<radius)
                 area=areatriangle(A,B,C);
             else
                 area=arealens(center, radius, B, A);
             end
        elseif (dB==radius)
             if (norm(C-center)<radius)
                 if (norm(list1(1,:)-B)<=norm(list1(2,:)-B))  
                    l11=list1(2,:);
                 else
                    l11=list1(1,:);                  
                  end
                  area=arealens(center, radius, list3(1,:), l11)+areatriangle(C,list3(1,:), l11)+areatriangle(B,C, l11);
             else
                 [b]=is_in_same_half_space(C,center,B, A);
                 if (b==1)
                     area=pi*(radius^2)-arealens(center,radius,list1(1,:),list1(2,:));
                 else
                     area=arealens(center,radius,list1(1,:),list1(2,:));
                 end
             end
        elseif (dA==radius)
            if (norm(C-center)<radius)
                 if (norm(list1(1,:)-A)<=norm(list1(2,:)-A))  
                    l11=list1(2,:);
                   else
                    l11=list1(1,:);                  
                  end
                  area=arealens(center, radius, list2(1,:), l11)+areatriangle(C,list2(1,:), l11)+areatriangle(C,A,l11);
             else
                 [b]=is_in_same_half_space(C,center,B,A);
                 if (b==1)
                     area=pi*(radius^2)-arealens(center,radius,list1(1,:),list1(2,:));
                 else
                     area=arealens(center,radius,list1(1,:),list1(2,:));
                 end
             end
        end         
    elseif (code == 17)
        %122
        if (norm(C-list2(1,:))<norm(C-list2(2,:)))
            l21=list2(1,:);
            l22=list2(2,:);
        else
            l21=list2(2,:);
            l22=list2(1,:);
        end
        
        if (norm(A-list1(1,:))<norm(A-list1(2,:)))
            l11=list1(1,:);
            l12=list1(2,:);
        else
            l11=list1(2,:);
            l12=list1(1,:);
        end
        
        if ((dA > radius) && (dB > radius) && (dC > radius))
            [b]=is_in_same_half_space(B,center,l11,l21);
            if (b==1)
               area =arealens(center, radius, l12, l22)+arealens(center, radius, l11, l21)+areatriangle(l12, l22, l11)+areatriangle(l22, l11, l21);
            else
               area =arealens(center, radius, l12, l22)+pi*(radius^2)-arealens(center, radius, l11, l21)+areatriangle(l12, l22, l11)+areatriangle(l22, l11, l21); 
            end
        elseif ((dA == radius) && (dB==radius))
            area =arealens(center, radius,A, l21)+areatriangle(B,A,l21);
        elseif ((dB == radius) && (dC==radius))
            area =arealens(center, radius,C, l11)+areatriangle(B,C,l11);
        elseif (dB==radius)
            [b]=is_in_same_half_space(B,center,l11,l21);
            if (b==1)
                area =arealens(center, radius,l21, l11)+areatriangle(B,l11,l21);
            else
                area =pi*(radius^2)-arealens(center, radius,l21, l11)+areatriangle(B,l11,l21);
            end
        elseif (dC==radius)
             area =arealens(center, radius,C,l11)+arealens(center, radius,l12,l22)+areatriangle(C,l11,l12)+areatriangle(C,l12,l22);
        elseif (dA==radius)
             area =arealens(center,radius,A,l21)+arealens(center, radius,l12,l22)+areatriangle(A,l21,l22)+areatriangle(A,l22,l12);
        end
    elseif (code == 23)
        %212
        if (norm(B-list1(1,:))<norm(B-list1(2,:)))
            l11=list1(1,:);
            l12=list1(2,:);
        else
            l11=list1(2,:);
            l12=list1(1,:);
        end
        
        if (norm(A-list3(1,:))<norm(A-list3(2,:)))
            l31=list3(1,:);
            l32=list3(2,:);
        else
            l31=list3(2,:);
            l32=list3(1,:);
        end
        
        if ((dA > radius) && (dB > radius) && (dC > radius))
            [b]=is_in_same_half_space(A,center,l11,l32);
            if (b==1)
                area =arealens(center, radius, l32, l11)+arealens(center, radius, l31, l12)+areatriangle(l11, l32, l12)+areatriangle(l32, l12, l31);
            else
                area =pi*(radius^2)-arealens(center, radius, l32, l11)+arealens(center, radius, l31, l12)+areatriangle(l11, l32, l12)+areatriangle(l32, l12, l31);
            end
        elseif ((dA == radius) && (dC==radius))
            area =arealens(center, radius,C, l11)+areatriangle(C,A,l11);
        elseif ((dA == radius) && (dB==radius))
            area =arealens(center, radius,B, l32)+areatriangle(B,A,l32);
        elseif (dA==radius)
            [b]=is_in_same_half_space(A,center,l32,l11);
            if (b==1)
                area =arealens(center, radius,l32, l11)+areatriangle(A,l11,l32);
            else
                area =pi*(radius^2)-arealens(center, radius,l32, l11)+areatriangle(A,l11,l32);
            end
        elseif (dB==radius)
             area =arealens(center, radius,B,l32)+arealens(center, radius,l31,l12)+areatriangle(B,l32,l31)+areatriangle(B,l31,l12);
        elseif (dC==radius)
             area =arealens(center,radius,C,l11)+arealens(center, radius,l31,l12)+areatriangle(C,l11,l31)+areatriangle(l11,l31,l12);
        end
    elseif (code == 25)
        %221
        if (norm(C-list2(1,:))<norm(C-list2(2,:)))
            l21=list2(1,:);
            l22=list2(2,:);
        else
            l21=list2(2,:);
            l22=list2(1,:);
        end
        
        if (norm(A-list3(1,:))<norm(A-list3(2,:)))
            l31=list3(1,:);
            l32=list3(2,:);
        else
            l31=list3(2,:);
            l32=list3(1,:);
        end
        
        if ((dA > radius) && (dB > radius) && (dC > radius))
            [b]=is_in_same_half_space(C,center,l31,l22);
            if (b==1)
                area =arealens(center, radius, l32, l21)+arealens(center, radius, l31, l22)+areatriangle(l32, l21, l31)+areatriangle(l21, l31, l22);
            else
                area =arealens(center, radius, l32, l21)+pi*(radius^2)-arealens(center, radius, l31, l22)+areatriangle(l32, l21, l31)+areatriangle(l21, l31, l22);
            end
        elseif ((dA == radius) && (dC==radius))
            area =arealens(center, radius,A, l22)+areatriangle(C,A,l22);
        elseif ((dB == radius) && (dC==radius))
            area =arealens(center, radius,B, l31)+areatriangle(B,C,l31);
        elseif (dC==radius)
            [b]=is_in_same_half_space(C,center,l22,l31);
            if (b==1)
                area =arealens(center, radius,l22, l31)+areatriangle(C,l22,l31);
            else
                area =pi*(radius^2)-arealens(center, radius,l22, l31)+areatriangle(C,l22,l31);
            end
         elseif (dB==radius)
             area =arealens(center, radius,B,l31)+arealens(center, radius,l32,l21)+areatriangle(B,l32,l31)+areatriangle(B,l32,l21);
         elseif (dA==radius)
             area =arealens(center,radius,A,l22)+arealens(center, radius,l32,l21)+areatriangle(A,l21,l22)+areatriangle(A,l21,l32);
        end   
    elseif (code == 26)
        %222
        
        if (norm(A-list1(1,:))<norm(A-list1(2,:)))
            l11=list1(1,:);
            l12=list1(2,:);
        else
            l11=list1(2,:);
            l12=list1(1,:);
        end
        
        if (norm(B-list2(1,:))<norm(B-list2(2,:)))
            l21=list2(1,:);
            l22=list2(2,:);
        else
            l21=list2(2,:);
            l22=list2(1,:);
        end
        
        if (norm(C-list3(1,:))<norm(C-list3(2,:)))
            l31=list3(1,:);
            l32=list3(2,:);
        else
            l31=list3(2,:);
            l32=list3(1,:);
        end
        
        if ((dA > radius) && (dB > radius) && (dC > radius))
            area = arealens(center, radius, l32, l11)+arealens(center, radius, l31, l22)+arealens(center, radius, l21, l12)+areatriangle(l32, l11, l31)+areatriangle(l11, l31, l22)+areatriangle(l11, l22, l12)+areatriangle(l12, l22, l21);
                
        elseif ((dA == radius) && (dB == radius) && (dC == radius))
            area=areatriangle(A,B,C);             
        elseif ((dA == radius) && (dC == radius) && (dB > radius))
            area =  arealens(center, radius, l12, l21)+areatriangle(A,C,l21)+areatriangle(A,l21,l12);
                
        elseif ((dC == radius) && (dB == radius) && (dA > radius))
            area =arealens(center, radius, l11, l32)+areatriangle(C,B,l32)+areatriangle(B,l11,l32);
        
        elseif ((dB == radius) && (dA == radius) && (dC > radius))
            area = arealens(center, radius, l22, l31)+areatriangle(B,A,l22)+areatriangle(A,l22,l31);
        elseif ((dA > radius) && (dC == radius) && (dB > radius))
            area =arealens(center, radius, l11, l32)+arealens(center, radius, l21, l12)+areatriangle(C,l21, l12)+areatriangle(C,l12, l11)+areatriangle(C,l11,l32);

        elseif ((dC > radius) && (dB == radius) && (dA > radius))
            area = arealens(center, radius, l22, l31)+arealens(center, radius, l32, l11)+areatriangle(B, l22, l31)+areatriangle(B, l31, l32)+areatriangle(B, l32, l11);
                
        elseif ((dB > radius) && (dA == radius) && (dC > radius))
            area =  arealens(center, radius,l12,l21)+arealens(center, radius,l31,l22)+areatriangle(A,l31, l22)+areatriangle(A,l22,l21)+areatriangle(A,l21,l12);
            
        end
        
    end
            
        
