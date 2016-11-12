clear all; close all; clc;

u_init = -3;
up_init = 3;

L = 10;
h = .125; %step size
x = -L:h:L;
y = x;

[X,Y] = meshgrid(x,y); 
U = zeros(161,161);
m = (length(x)-1)/4;
x1 = linspace(-L,L,m);
y1 = x1;
[X1,Y1] = meshgrid(x1,y1);      %x and y are small grid, indexed from 1 to 161
                                %x1 and y1 and big grid, indexed from 1 to
                                %40
                               
%step 1
point = [0,0];
goal_ind = [randi(m) , randi(m)];
goal = [x1(goal_ind(1)),y1(goal_ind(2))];
remove = zeros(m^2,2);
remove(1,:) = goal_ind;
remove_comp = zeros(m^2,2);
for i = 1:m
    for j = 1:m
        remove_comp(m*i+j,1) = i;
        remove_comp(m*i+j,2) = j;
    end
end

h = .3; %resets step size to what works well for paths on the grid

[ap, bp] = pades(point(1)+point(2)*1i, u_init, up_init); %find pade coeffs at ORIGIN
points_dir = cell(1,5);
points_dir{1}(1,:) = point;
points_dir{2}(1,1) = u_init;
points_dir{3}(1,1) = up_init;
points_dir{4}(1,:) = ap;
points_dir{5}(1,:) = bp;
countp = 2;
old_countp = 1;         % these are for
point1 = [0,0];         %plotting later on


for k=1:m*m;

    while norm(point-goal) >= h

       vec = (goal-point);
       vec = vec/norm(vec);
       vec = h*vec;
       options = zeros(5,2); %listing of step options
       options(1,:) = vec;
       theta = 22.5;
       R1 = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
       options(2,:) = vec*R1;
       theta2 = -22.5;
       R2 = [cosd(theta2) -sind(theta2); sind(theta2) cosd(theta2)];
       options(3,:) = vec*R2;
       theta3 = 45;
       R3 = [cosd(theta3) -sind(theta3); sind(theta3) cosd(theta3)];
       options(4,:) = vec*R3;
       theta4 = -45;
       R4 = [cosd(theta4) -sind(theta4); sind(theta4) cosd(theta4)];
       options(5,:) = vec*R4;
       
       options1 = bsxfun(@plus,options,point); %this sets the possible points FROM 
                                               % the original point!!!
       % the following 2 lines evaluate and choose the minimum abs
       cases = polyval(ap, options(:,1) + 1i*options(:,2))...
           ./polyval(bp, options(:,1) + 1i*options(:,2));
       [ab, ind] = min( abs( cases ) );

       u = cases(ind);       %sets u for correct point
       point = [options1(ind,1) , options1(ind,2)];
       ht = options(ind,1) + 1i*options(ind,2);
       z = point(1) + 1i*point(2);   %sets z value for best point
       [app, bpp] = polyder(ap,bp);             %derivative of pade
       up = polyval(app,ht)/polyval(bpp,ht);    %finds u' for correct point
       [ap, bp] = pades(z, u, up); %find pade coeffs at POINT

       %following lines store the selected point into the cell structure
       points_dir{1}(countp,:) = point;
       points_dir{2}(countp,1) = u;
       points_dir{3}(countp,1) = up;
       points_dir{4}(countp,:) = ap;
       points_dir{5}(countp,:) = bp;
       countp = countp+1;
    end
    % following plot command produces pretty picture! old_countp is for
    % different branches, point1 is to include the previous point that is
    % searched for
     plot([point1(1,1); points_dir{1}(old_countp:countp-1 ,1)],[point1(1,2); points_dir{1}(old_countp:countp-1, 2)]);
     hold on;
     old_countp = countp;
    
     v = setdiff(remove_comp, remove,'rows');                 %these lines remove the
     goal_ind(1) = v(randi(numel(v(:,1) )));            %achieved goal and set a new one
     goal_ind(2) = v(randi(numel(v(:,2) ))); 
     
     goal = [x1(goal_ind(1)),y1(goal_ind(2))];
     remove(k+1,:) = goal_ind;
     
     %find the closest point (bsxfun does row-wise subtraction)
     [dist, indd]= min( sqrt(sum(abs( bsxfun(@minus,goal,points_dir{1}) ).^2,2)) );
     point = points_dir{1}(indd,:);
     u = points_dir{2}(indd,1);
     up = points_dir{3}(indd,1);
     ap = points_dir{4}(indd,:);
     bp = points_dir{5}(indd,:);
     point1 = point;                %for plot purposes

end

%This marks the end of "step one." 
%Step 2 is to use all of this information to find the points on the fine
%grid
for i = 1:length(x)
    for j = 1:length(y)
    [dist, indg]= min(sqrt(sum(abs( bsxfun(@minus, [x(i),y(j)],points_dir{1})).^2,2)));
    z = x(i)+1i*y(j);
    ht = z - points_dir{1}(indg,1) + points_dir{1}(indg,2)*1i;
    U(i,j) = polyval(points_dir{4}(indg,:),ht)/ polyval(points_dir{5}(indg,:),ht);
        
    end
end






