function out=brick8(mode,b,c,d,e)
  
% BEAM3 does as listed below. It is an Euler-Bernoulli
% beam/rod/torsion model. 
% Beam properties (bprops) are in the order
% bprops=[E G rho A1 A2 A3 J1 J2 J3 Ixx1 Ixx2 Ixx3 Iyy1 Iyy2 Iyy3]
% Third node is in the middle.
% Fourth "node" defines the beam y plane and is actually from the
% points array.
%%
% Defining beam element properties in wfem input file:
% element properties
%   E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3 
% Torsional rigidity, $J$, must be less than or equal
% to $Iyy+Izz$ at any given cross section.  
%
% Defining beam3 element in wfem input file:
%   node1 node2 node3 pointnumber materialnumber 

%
% See wfem.m for more explanation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Variables (global):
% -------------------
% K       :    Global stiffness matrix
% Ks      :    Global stiffness buckling matrix
% M       :    Global mass matrix
% nodes   :    [x y z] nodal locations
global ismatnewer
global K
global Ks
global M
global nodes % Node locations
global elprops
global element
global points
global Fepsn % Initial strain "forces". 
global lines
global restart
global reload
global curlineno
global DoverL
global surfs
%
% Variables (local):
% ------------------
% bnodes  :    node/point numbers for actual beam nodes 1-2-3 and point
% k       :    stiffness matrix in local coordiates
% kg      :    stiffness matrix rotated into global coordinates
% m       :    mass matrix in local coordiates
% mg      :    mass matrix rotated into global coordinates
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright Joseph C. Slater, 7/26/2002.
% joseph.slater@wright.edu
out=0;
if strcmp(mode,'numofnodes')
    % This allows a code to find out how many nodes this element has
    out=8;
end
if strcmp(mode,'generate')
  elnum=c;%When this mode is called, the element number is the 3rd
          %argument. 
  
          %The second argument (b) is the element
          %definition. For this element b is
          %node1 node2 node3 node4 node5 node6 node7 node8 and material#
  
          %There have to be 9 numbers for this element's
          %definition (above)
  if length(b)==9
      element(elnum).nodes=b(1:8);
%       %If the user puts the middle node in the wrong place, tell them.
%       if norm((nodes(b(1),:)+nodes(b(2),:))/2-nodes(b(3),:))/ ...
%               norm(nodes(b(1),:)-nodes(b(2),:))>.001
%           disp(['WARNING: Node ' num2str(b(3)) ...
%                 ' is not in the middle of element ' ...
%                 num2str(elnum) ' on line ' num2str(curlineno) '.'])
%       end
      element(elnum).properties=b(9);
%       element(elnum).point=b(4);
  else 
	  b
      %There have to be 9 numbers on a line defining the
      %element. 
      warndlg(['Element ' num2str(elnum) ' on line ' ...
               num2str(element(elnum).lineno) ' entered incorrectly.'], ...
              ['Malformed Element'],'modal')
      return
  end
 
end

% Here we figure out what the beam properties mean. If you need
% them in a mode, that mode should be in the if on the next line.
if strcmp(mode,'make')||strcmp(mode,'istrainforces')
  elnum=b;% When this mode is called, the element number is given
          % as the second input.
  bnodes=[element(elnum).nodes];% The point is
                                                     % referred to
                                                     % as node 4
                                                     % below,
                                                     % although it
                                                     % actually
                                                     % calls the
                                                     % array points
                                                     % to get its
                                                     % location. Its
                                                     % not really a
                                                     % node, but
                                                     % just a point
                                                     % that helps
                                                     % define
                                                     % orientation. Your
                                                     % element may
                                                     % not need
                                                     % such a
                                                     % reference point.
  bprops=elprops(element(elnum).properties).a;% element(elnum).properties 
                                              % stores the
                                              % properties number
                                              % of the current
                                              % elnum. elprops
                                              % contains this
                                              % data. This is
                                              % precisely the
                                              % material properties
                                              % line in an
                                              % array. You can pull
                                              % out any value you
                                              % need for your use. 
  
% 
  if length(bprops)==3
      E=bprops(1);
      G=bprops(2);
      rho=bprops(3);
%       A1=bprops(4);
%       A2=bprops(5);
%       A3=bprops(6);
%       J1=bprops(7);
%       J2=bprops(8);
%       J3=bprops(9);
%       Izz1=bprops(10);
%       Izz2=bprops(11);
%       Izz3=bprops(12);
%       Iyy1=bprops(13);
%       Iyy2=bprops(14);
%       Iyy3=bprops(15);
  else
      warndlg(['The number of material properties set for ' ...
               'this element (' num2str(length(bprops)) ') isn''t ' ...
               'appropriate for a beam3 element. '      ...
               'Please refer to the manual.'],...
              'Bad element property definition.','modal');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beam properties (bprops) are in the order
% bprops=[E G rho A1 A2 A3 J1 J2 J3 Izz1 Izz2 Izz3 Iyy1 Iyy2 Iyy3]
% For a linear beam they are
% bprops=[E G rho A1 A2 J1 J2 Izz1 Izz2 Iyy1 Iyy2]

if strcmp(mode,'make')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Define beam node locations for easy later referencing
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  x1=nodes(bnodes(1),1);
  y1=nodes(bnodes(1),2);
  z1=nodes(bnodes(1),3);
  x2=nodes(bnodes(2),1);
  y2=nodes(bnodes(2),2);
  z2=nodes(bnodes(2),3);
  x3=nodes(bnodes(3),1);
  y3=nodes(bnodes(3),2);
  z3=nodes(bnodes(3),3);
  x4=nodes(bnodes(4),1);
  y4=nodes(bnodes(4),2);
  z4=nodes(bnodes(4),3);
  x5=nodes(bnodes(5),1);
  y5=nodes(bnodes(5),2);
  z5=nodes(bnodes(5),3);
  x6=nodes(bnodes(6),1);
  y6=nodes(bnodes(6),2);
  z6=nodes(bnodes(6),3);
  x7=nodes(bnodes(7),1);
  y7=nodes(bnodes(7),2);
  z7=nodes(bnodes(7),3);
  x8=nodes(bnodes(8),1);
  y8=nodes(bnodes(8),2);
  z8=nodes(bnodes(8),3);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Shape functions for higher order beam. 
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Shape functions in matrix polynomial form (polyval style) for brick8
  Results = brick8_ShapeFun();
  bn{1} = Results.c1;
  bnd{1} = Results.c1d;
  
  bn{2} = Results.c2;
  bnd{2} = Results.c2d; 
  
  bn{3} = Results.c3;
  bnd{3} = Results.c3d;

  bn{4} = Results.c4;
  bnd{4} = Results.c4d;  
  
  bn{5} = Results.c5;
  bnd{5} = Results.c5d;
  
  bn{6} = Results.c6;
  bnd{6} = Results.c6d;
  
  bn{7} = Results.c7;
  bnd{7} = Results.c7d;
  
  bn{8} = Results.c8;
  bnd{8} = Results.c8d;  
  
  
  numbeamgauss=7; % Number of Gauss points for integration in the 1D case (such as a beam)
  [bgpts,bgpw]=gauss(numbeamgauss);
  kb=zeros(24,24);% For this brick8 element, 8 nodes, 3DOF each, is a 24 by 24
                  % matrix. 


  X = [x1, x2, x3, x4, x5, x6, x7, x8];
  Y = [y1, y2, y3, y4, y5, y6, y7, y8];
  Z = [z1, z2, z3, z4, z5, z6, z7, z8];


  for i = 1:numbeamgauss
      for j = 1:numbeamgauss
          for k = 1:numbeamgauss
              % Coordinate of Gauss Point in a brick in the form [x,y,z]
              gpts = [bgpts(i),bgpts(j),bgpts(k)];  
              % Evaluate the Jocobian at Gauss Point
              J = jacobian(X,Y,Z,gpts); 
              J_det = det(J);
              B = [];
              % Calculate dN_i/dX, dN_i/dY, dN_i/dZ
              for n=1:8
                  % isoparametric coordinate
                  dNdx = poly3dval(bnd{n}(1,:),gpts);
                  dNdy = poly3dval(bnd{n}(2,:),gpts);
                  dNdz = poly3dval(bnd{n}(3,:),gpts);
                  % Cartesion Coordinate
                  dNdXYZ = J\[dNdx;dNdy;dNdz];
                  Bi = [dNdXYZ(1), 0, 0
                          0, dNdXYZ(2), 0
                          0, 0, dNdXYZ(3)
                          dNdXYZ(2), dNdXYZ(1), 0
                          0, dNdXYZ(3), dNdXYZ(2)
                          dNdXYZ(3), 0, dNdXYZ(1)];
                  B = [B,Bi];
              end  
              kb = kb + bgpw(i)*bgpw(j)*bgpw(k)*B'*E*B*J_det;
          end
      end
  end
  [vectors,values]=eig(kb);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  % Derivation of Mass matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  numbeamgauss=numbeamgauss+3; %Need more gauss points for the mass
                               %matrix. 
  [bgpts,bgpw]=gauss(numbeamgauss);
  mb=zeros(24,24); %initialize empty mass matrix
  
  for i = 1:numbeamgauss
      for j = 1:numbeamgauss
          for k = 1:numbeamgauss
              % Coordinate of Gauss Point in a brick in the form [x,y,z]
              gpts = [bgpts(i),bgpts(j),bgpts(k)];  
              % Evaluate the Jocobian at Gauss Point
              J = jacobian(X,Y,Z,gpts); 
              J_det = det(J);
              N = [];
              for n = 1:8
                  Ni = poly3dval(bn{n},gpts);
                  N = [N, diag([Ni, Ni, Ni])];
              end                       
              mb = mb + bgpw(i)*bgpw(j)*bgpw(k)*N'*rho*N*J_det;
          end
      end
  end
  mb;
   
     
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Coordinate rotations
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  R1=([x2 y2 z2]-[x1 y1 z1]);% Vector along element
  lam1=R1/norm(R1);% Unit direction
  R2=([x4 y4 z4]-[x1 y1 z1]);%  Vector to the point
  R2perp=R2-dot(R2,lam1)*lam1;% Part of R2 perpendicular to lam1
  udirec=0;
  while norm(R2perp)<10*eps% If R2perp is too small, (point in line
                           % with element, we need to cover the
                           % users a$$ and generate a point that
                           % isn't. We should put out a warning,
                           % but I commented it out. 
    udirec=udirec+1;
    disp('oops, point is on the line of the element'); %This was my warning. 
    %pause
    [minval,minloc]=min(lam1);
    R2perp=zeros(1,3);
    R2perp(udirec)=1;
    R2perp=R2perp-dot(R2perp,lam1)*lam1;
  end
  %Make the unit direction vectors for rotating and put them in the
  %rotation matrix. 
  lam2=R2perp/norm(R2perp);
  lam3=cross(lam1,lam2);
  lamloc=[lam1;lam2;lam3];
  lam=sparse(18,18);
  lam(1:3,1:3)=lamloc;
  lam(4:6,4:6)=lamloc;
  lam(7:9,7:9)=lamloc;
  lam(10:12,10:12)=lamloc;
  lam(13:15,13:15)=lamloc;
  lam(16:18,16:18)=lamloc;
  
% $$$     lam=[lamloc z z z z z;
% $$$          z lamloc z z z z;
% $$$          z z lamloc z z z;
% $$$          z z z lamloc z z;
% $$$          z z z z lamloc z;
% $$$          z z z z z lamloc];
  element(elnum).lambda=lam;
  element(elnum).m=m;
  element(elnum).k=k;

  kg=lam'*k*lam;
  mg=lam'*m*lam;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Assembling matrices into global matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  bn1=bnodes(1);bn2=bnodes(2);bn3=bnodes(3);
  indices=[bn1*6+(-5:0) bn2*6+(-5:0) bn3*6+(-5:0)] ;

  K(indices,indices)=K(indices,indices)+kg;
  M(indices,indices)=M(indices,indices)+mg;

  % At this point we also know how to draw the element (what lines
  % and surfaces exist). For the beam3 element, 2 lines are
  % appropriate. Just add the pair of node numbers to the lines
  % array and that line will always be drawn.
  numlines=size(lines,1);
  lines(numlines+1,:)=[bn1 bn3];
  lines(numlines+2,:)=[bn3 bn2];
  
  % If I have 4 nodes that I want to use to represent a surface, I
  % do the following.
  panelcolor=[1 0 1];% This picks a color. You can change the
                     % numbers between 0 and 1. 
  %Don't like this color? Use colorui to pick another one. Another
  %option is that if we can't see the elements separately we can
  %chunk up x*y*z, divide by x*y*x of element, see if we get
  %integer powers or not to define colors that vary by panel. 
  
  
  % You need to uncomment this line and assign values to node1,
  % node2, node3, and node4 in order to draw A SINGLE SURFACE. For
  % a brick, you need 6 lines like this. 
  %surfs=[surfs;node1 node2 node3 node4 panelcolor];
  
  %Each surface can have a different color if you like. Just change
  %the last three numbers on the row corresponding to that
  %surface. 

elseif strcmp(mode,'istrainforces')
  % You don't need this
  % We need to have the stiffness matrix and the coordinate roation matrix.
 
elseif strcmp(mode,'draw')
elseif strcmp(mode,'buckle')
end