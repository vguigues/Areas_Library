function triangles = polygon_triangulate ( n, x, y )

%*****************************************************************************80
%
%% POLYGON_TRIANGULATE determines a triangulation of a polygon.
%
%  Discussion:
%
%    There are N-3 triangles in the triangulation.
%
%    For the first N-2 triangles, the first edge listed is always an
%    internal diagonal.
%
%    Thanks to Gene Dial for pointing out a mistake in the area calculation,
%    10 September 2016.
%-
%    Gene Dial requested an angle tolerance of about 1 millionth radian or 
%    5.7E-05 degrees, 26 June 2018.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    26 June 2018
%
%  Author:
%
%    Original C version by Joseph ORourke.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Joseph ORourke,
%    Computational Geometry in C,
%    Cambridge, 1998,
%    ISBN: 0521649765,
%    LC: QA448.D38.
%
%  Parameters:
%
%    Input, integer N, the number of vertices.
%
%    Input, real X(N), Y(N), the coordinates of each vertex.
%
%    Output, integer TRIANGLES(3,N-2), the triangles of the triangulation.
%
  angle_tol = 5.7E-05;
%
%  We must have at least 3 vertices.
%
  if ( n < 3 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'POLYGON_TRIANGULATE - Fatal error!\n' );
    fprintf ( 1, '  N < 3.\n' );
    error ( 'POLYGON_TRIANGULATE - Fatal error!' );
  end
%
%  Consecutive vertices cannot be equal.
%
  node_m1 = n;
  for node = 1 : n
    if ( x(node_m1) == x(node) && y(node_m1) == y(node) )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'POLYGON_TRIANGULATE - Fatal error!\n' );
      fprintf ( 1, '  Two consecutive nodes are identical.\n' );
      error ( 'POLYGON_TRIANGULATE - Fatal error!' );
    end
    node_m1 = node;
  end
%
%  No node can be the vertex of an angle less than 1 degree 
%  in absolute value.
%
  node1 = n;

  for node2 = 1 : n

    node3 = mod ( node2, n ) + 1;

    angle = angle_degree ( ...
      x(node1), y(node1), ...
      x(node2), y(node2), ...
      x(node3), y(node3) );

    if ( abs ( angle ) <= angle_tol )
      fprintf ( 1, '\n' );
      fprintf ( 1, 'POLYGON_TRIANGULATE - Fatal error!\n' );
      fprintf ( 1, '  Polygon has an angle smaller than %g.\n', angle_tol );
      fprintf ( 1, '  occurring at node %d\n', node2 );
      error ( 'POLYGON_TRIANGULATE - Fatal error!' );
    end

    node1 = node2;

  end
%
%  Area must be positive.
%
  area = polygon_area ( n, x, y );

  if ( area <= 0.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'POLYGON_TRIANGULATE - Fatal error!\n' );
    fprintf ( 1, '  Polygon has zero or negative area.\n' );
    error ( 'POLYGON_TRIANGULATE - Fatal error!' );
  end
%
%  PREV_NODE and NEXT_NODE point to the previous and next nodes.
%
  prev_node = zeros ( n, 1 );
  next_node = zeros ( n, 1 );

  i = 1;
  prev_node(i) = n;
  next_node(i) = i + 1;

  for i = 2 : n - 1
    prev_node(i) = i - 1;
    next_node(i) = i + 1;
  end

  i = n;
  prev_node(i) = i - 1;
  next_node(i) = 1;
%
%  EAR indicates whether the node and its immediate neighbors form an ear
%  that can be sliced off immediately.
%
  ear = zeros ( n, 1 );

  for i = 1 : n
    ear(i) = diagonal ( prev_node(i), next_node(i), n, prev_node, ...
      next_node, x, y );
  end

  first = 1;

  triangles = zeros ( 3, n - 2 );
  triangle_num = 0;

  i2 = first;

  while ( triangle_num < n - 3 )
%
%  If I2 is an ear, gather information necessary to carry out
%  the slicing operation and subsequent "healing".
%
    if ( ear(i2) )

      i3 = next_node(i2);
      i4 = next_node(i3);
      i1 = prev_node(i2);
      i0 = prev_node(i1);
%
%  Make vertex I2 disappear.
%
      next_node(i1) = i3;
      prev_node(i3) = i1;
%
%  Update the earity of I1 and I3, because I2 disappeared.
%
      ear(i1) = diagonal ( i0, i3, n, prev_node, next_node, x, y );
      ear(i3) = diagonal ( i1, i4, n, prev_node, next_node, x, y );
%
%  Add the diagonal [I3, I1, I2] to the list.
%
      triangle_num = triangle_num + 1;
      triangles(1,triangle_num) = i3;
      triangles(2,triangle_num) = i1;
      triangles(3,triangle_num) = i2;

    end
%
%  Try the next vertex.
%
    i2 = next_node(i2);

  end
%
%  The last triangle is formed from the three remaining vertices.
%
  i3 = next_node(i2);
  i1 = prev_node(i2);

  triangle_num = triangle_num + 1;
  triangles(1,triangle_num) = i3;
  triangles(2,triangle_num) = i1;
  triangles(3,triangle_num) = i2;

  return
end
function value = angle_degree ( x1, y1, x2, y2, x3, y3 )

%*****************************************************************************80
%
%% ANGLE_DEGREE returns the degree angle defined by three points.
%
%  Discussion:
%
%        P1
%        /
%       /
%      /
%     /
%    P2--------->P3
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    28 August 2016
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real X1, Y1, X2, Y2, X3, Y3, the coordinates of the points
%    P1, P2, P3.
%
%    Output, real VALUE, the angle swept out by the rays, measured
%    in degrees.  0 <= VALUE < 360.  If either ray has zero length,
%    then VALUE is set to 0.
%
  x = ( x3 - x2 ) * ( x1 - x2 ) + ( y3 - y2 ) * ( y1 - y2 );

  y = ( x3 - x2 ) * ( y1 - y2 ) - ( y3 - y2 ) * ( x1 - x2 );

  if ( x == 0.0 && y == 0.0 )
    value = 0.0;
    return
  end

  value = atan2 ( y, x );

  if ( value < 0.0 )
    value = value + 2.0 * pi;
  end

  value = 180.0 * value / pi;

  return
end
function value = between ( xa, ya, xb, yb, xc, yc )

%*****************************************************************************80
%
%% BETWEEN is TRUE if vertex C is between vertices A and B.
%
%  Discussion:
%
%    The points must be (numerically) collinear.
%
%    Given that condition, we take the greater of XA - XB and YA - YB
%    as a "scale" and check where C's value lies.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    03 May 2014
%
%  Author:
%
%    Original C version by Joseph ORourke.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Joseph ORourke,
%    Computational Geometry in C,
%    Cambridge, 1998,
%    ISBN: 0521649765,
%    LC: QA448.D38.
%
%  Parameters:
%
%    Input, real XA, YA, XB, YB, XC, YC, the coordinates of the vertices.
%
%    Output, integer VALUE, is 1 if C is between A and B,
%    and 0 otherwise.
%
  if ( ~ collinear ( xa, ya, xb, yb, xc, yc ) )
    value = 0;
  elseif ( abs ( ya - yb ) < abs ( xa - xb ) )
    xmax = max ( xa, xb );
    xmin = min ( xa, xb );
    value = ( xmin <= xc && xc <= xmax );
  else
    ymax = max ( ya, yb );
    ymin = min ( ya, yb );
    value = ( ymin <= yc && yc <= ymax );
  end

  return
end
function value = collinear ( xa, ya, xb, yb, xc, yc )

%*****************************************************************************80
%
%% COLLINEAR returns a measure of collinearity for three points.
%
%  Discussion:
%
%    In order to deal with collinear points whose coordinates are not
%    numerically exact, we compare the area of the largest square
%    that can be created by the line segment between two of the points
%    to (twice) the area of the triangle formed by the points.
%
%    If the points are collinear, their triangle has zero area.
%    If the points are close to collinear, then the area of this triangle
%    will be small relative to the square of the longest segment.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    10 September 2016
%
%  Author:
%
%    Original C version by Joseph ORourke.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Joseph ORourke,
%    Computational Geometry in C,
%    Cambridge, 1998,
%    ISBN: 0521649765,
%    LC: QA448.D38.
%
%  Parameters:
%
%    Input, real XA, YA, XB, YB, XC, YC, the vertex coordinates.
%
%    Output, integer VALUE, is 1 if the points are judged to be collinear,
%    and 0 otherwise.
%
  area = triangle_area ( xa, ya, xb, yb, xc, yc );

  side_ab_sq = ( xa - xb ) ^ 2 + ( ya - yb ) ^ 2;
  side_bc_sq = ( xb - xc ) ^ 2 + ( yb - yc ) ^ 2;
  side_ca_sq = ( xc - xa ) ^ 2 + ( yc - ya ) ^ 2;

  side_max_sq = max ( side_ab_sq, max ( side_bc_sq, side_ca_sq ) );

  if ( side_max_sq <= eps )
    value = 1;
  elseif ( 2.0 * abs ( area ) <= eps * side_max_sq )
    value = 1;
  else
    value = 0;
  end

  return
end
function value = diagonal ( im1, ip1, n, prev_node, next_node, x, y )

%*****************************************************************************80
%
%% DIAGONAL: VERTEX(IM1) to VERTEX(IP1) is a proper internal diagonal.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    03 May 2014
%
%  Author:
%
%    Original C version by Joseph ORourke.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Joseph ORourke,
%    Computational Geometry in C,
%    Cambridge, 1998,
%    ISBN: 0521649765,
%    LC: QA448.D38.
%
%  Parameters:
%
%    Input, integer IM1, IP1, the indices of two vertices.
%
%    Input, integer N, the number of vertices.
%
%    Input, integer PREV_NODE(N), the previous neighbor of each vertex.
%
%    Input, integer NEXT_NODE(N), the next neighbor of each vertex.
%
%    Input, real X(N), Y(N), the coordinates of each vertex.
%
%    Output, integer VALUE, the value of the test.
%
  value1 = in_cone ( im1, ip1, n, prev_node, next_node, x, y );
  value2 = in_cone ( ip1, im1, n, prev_node, next_node, x, y );
  value3 = diagonalie ( im1, ip1, n, next_node, x, y );

  value = value1 && value2 && value3;

  return
end
function value = diagonalie ( im1, ip1, n, next_node, x, y )

%*****************************************************************************80
%
%% DIAGONALIE: VERTEX(IM1):VERTEX(IP1) is a proper diagonal.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    08 February 2011
%
%  Author:
%
%    Original C version by Joseph ORourke.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Joseph ORourke,
%    Computational Geometry in C,
%    Cambridge, 1998,
%    ISBN: 0521649765,
%    LC: QA448.D38.
%
%  Parameters:
%
%    Input, integer IM1, IP1, the indices of two vertices.
%
%    Input, integer N, the number of vertices.
%
%    Input, integer NEXT_NODE(N), the next neighbor of each vertex.
%
%    Input, real X(N), Y(N), the coordinates of each vertex.
%
%    Output, int VALUE, the value of the test.
%
  first = im1;
  j = first;
  jp1 = next_node(first);

  value = 1;
%
%  For each edge VERTEX(J):VERTEX(JP1) of the polygon:

  while ( 1 )
%
%  Skip any edge that includes vertex IM1 or IP1.
%
    if ( j == im1 || j == ip1 || jp1 == im1 || jp1 == ip1 )

    else

      value2 = intersect ( x(im1), y(im1), x(ip1), y(ip1), x(j), y(j), ...
        x(jp1), y(jp1) );

      if ( value2 )
        value = 0;
        return
      end
    end

    j = jp1;
    jp1 = next_node(j);

    if ( j == first )
      break
    end

  end

  return
end
function value = in_cone ( im1, ip1, n, prev_node, next_node, x, y )

%*****************************************************************************80
%
%% IN_CONE is TRUE if the diagonal VERTEX(IM1):VERTEX(IP1) is strictly internal.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    06 February 2011
%
%  Author:
%
%    Original C version by Joseph ORourke.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Joseph ORourke,
%    Computational Geometry in C,
%    Cambridge, 1998,
%    ISBN: 0521649765,
%    LC: QA448.D38.
%
%  Parameters:
%
%    Input, integer IM1, IP1, the indices of two vertices.
%
%    Input, integer N, the number of vertices.
%
%    Input, integer PREV_NODE(N), the previous neighbor of each vertex.
%
%    Input, integer NEXT_NODE(N), the next neighbor of each vertex.
%
%    Input, real X(N), Y(N), the coordinates of each vertex.
%
%    Output, int VALUE, the value of the test.
%
  im2 = prev_node(im1);
  i = next_node(im1);

  t1 = triangle_area ( x(im1), y(im1), x(i), y(i), x(im2), y(im2) );

  if ( 0.0 <= t1 )

    t2 = triangle_area ( x(im1), y(im1), x(ip1), y(ip1), x(im2), y(im2) );
    t3 = triangle_area ( x(ip1), y(ip1), x(im1), y(im1), x(i), y(i) );

    value = ( ( 0.0 < t2 ) && ( 0.0 < t3 ) );

  else

    t4 = triangle_area ( x(im1), y(im1), x(ip1), y(ip1), x(i), y(i) );
    t5 = triangle_area ( x(ip1), y(ip1), x(im1), y(im1), x(im2), y(im2) );

    value = ~ ( ( 0.0 <= t4 ) && ( 0.0 <= t5 ) );

  end

  return
end
function value = intersect ( xa, ya, xb, yb, xc, yc, xd, yd )

%*****************************************************************************80
%
%% INTERSECT is true if lines VA:VB and VC:VD intersect.
%
%  Discussion:
%
%    Thanks to Gene Dial for correcting the call to intersect_prop,
%    08 September 2016.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    08 September 2016
%
%  Author:
%
%    Original C version by Joseph ORourke.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Joseph ORourke,
%    Computational Geometry in C,
%    Cambridge, 1998,
%    ISBN: 0521649765,
%    LC: QA448.D38.
%
%  Parameters:
%
%    Input, real XA, YA, XB, YB, XC, YC, XD, YD, the X and Y coordinates
%    of the four vertices.
%
%    Output, integer VALUE, the value of the test.
%
  if ( intersect_prop ( xa, ya, xb, yb, xc, yc, xd, yd ) )
    value = 1;
  elseif ( between ( xa, ya, xb, yb, xc, yc ) )
    value = 1;
  elseif ( between ( xa, ya, xb, yb, xd, yd ) )
    value = 1;
  elseif ( between ( xc, yc, xd, yd, xa, ya ) )
    value = 1;
  elseif ( between ( xc, yc, xd, yd, xb, yb ) )
    value = 1;
  else
    value = 0;
  end

  return
end
function value = intersect_prop ( xa, ya, xb, yb, xc, yc, xd, yd )

%*****************************************************************************80
%
%% INTERSECT_PROP is TRUE if lines VA:VB and VC:VD have a proper intersection.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    07 February 2011
%
%  Author:
%
%    Original C version by Joseph ORourke.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Joseph ORourke,
%    Computational Geometry in C,
%    Cambridge, 1998,
%    ISBN: 0521649765,
%    LC: QA448.D38.
%
%  Parameters:
%
%    Input, real XA, YA, XB, YB, XC, YC, XD, YD, the X and Y coordinates
%    of the four vertices.
%
%    Output, integer VALUE, the result of the test.
%
  if ( collinear ( xa, ya, xb, yb, xc, yc ) )
    value = 0;
  elseif ( collinear ( xa, ya, xb, yb, xd, yd ) )
    value = 0;
  elseif ( collinear ( xc, yc, xd, yd, xa, ya ) )
    value = 0;
  elseif ( collinear ( xc, yc, xd, yd, xb, yb ) )
    value = 0;
  else
    t1 = triangle_area ( xa, ya, xb, yb, xc, yc );
    t2 = triangle_area ( xa, ya, xb, yb, xd, yd );
    t3 = triangle_area ( xc, yc, xd, yd, xa, ya );
    t4 = triangle_area ( xc, yc, xd, yd, xb, yb );
    value = xor ( 0.0 < t1, 0.0 < t2 ) && ...
            xor ( 0.0 < t3, 0.0 < t4 );
  end

  return
end
function area = polygon_area ( n, x, y )

%*****************************************************************************80
%
%% POLYGON_AREA returns the area of a polygon.
%
%  Discussion:
%
%    The vertices should be listed in counter-clockwise order so that
%    the area will be positive.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    10 September 2016
%
%  Author:
%
%    John Burkardt.
%
%  Parameters:
%
%    Input, integer N, the number of vertices.
%
%    Input, real X(N), Y(N), the vertex coordinates.
%
%    Output, real AREA, the area of the polygon.
%
  area = 0.0;
  im1 = n;

  for i = 1 : n
    area = area + x(im1) * y(i) - x(i) * y(im1);
    im1 = i;
  end

  area = 0.5 * area;

  return
end
function area = triangle_area ( xa, ya, xb, yb, xc, yc )

%*****************************************************************************80
%
%% TRIANGLE_AREA computes the signed area of a triangle.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    07 February 2011
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real XA, YA, XB, YB, XC, YC, the vertices.
%
%    Output, real AREA, the signed area of the triangle.
%
  area = 0.5 * ( ( xb - xa ) * ( yc - ya ) ...
               - ( xc - xa ) * ( yb - ya ) );

  return
end
