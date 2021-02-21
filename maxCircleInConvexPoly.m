function [maxRad,maxRadOrigin]=maxCircleInConvexPoly(xPoly,yPoly,shrinkage)

% 07-17-20 Mitchell Tillman

%% WORKS FOR 2D CONVEX SHAPES ONLY!!!
%
%% Overview:
% "Inflates" (translates and increases  size) multiple "balloons" (circles) through an iterative process
% incrementing along the vectors that bisect two sides until largest possible radius is reached.
%
%% Inputs:
% xPoly: Vector of x-coordinates.
% yPoly: Vector of y-coordinates.
% Shrinkage: Argument to Matlab's 'boundary' function. 0 for maximal size hull.
%
%% Outputs:
% maxRad=radius of largest possible circle.
% maxRadOrigin=(x,y) origin of largest possible circle in convex polygon.

%% Specify plotting preferences.
plotRotOn=0; % Change to 1 to plot rotated polygon, otherwise 0.
plotOrigOn=0; % Change to 1 to plot original polygon, otherwise 0.
plotPrevBals=0; % Change to 1 to plot circles prior to largest, otherwise 0. Must also have plotRotOn or plotOrigOn set to 1

%% 1. Compute the boundary of the points.
if size(xPoly,1)<size(xPoly,2)
    xPoly=xPoly'; % Make column vector
end
if size(yPoly,1)<size(yPoly,2)
    yPoly=yPoly'; % Make column vector
end
poly=[xPoly yPoly]; % Aggregate X & Y data into columns of matrix.
[k,~]=boundary(poly,shrinkage); % Find order of points around polygon of given shrinkage.
boundPoly=poly(k,:); % Re-order the polygon points.

if plotOrigOn==1
    Q=figure;
    plot(boundPoly(:,1),boundPoly(:,2),'k');
    hold on;
    axis equal;
end

numSides=size(boundPoly,1)-1;

%% 2. Find the interior angle between each of the sides.
% Where sideI and sideII are both vectors.
% theta=acos(dot(sideI,sideII)/(norm(sideI)*norm(sideII)));
sideLengths=zeros(numSides,1); angleThetas=zeros(numSides,1);
for sideNum=1:numSides % Iterate through each pair of sides.
    if sideNum<numSides
        sideI=boundPoly(sideNum+1,:)-boundPoly(sideNum,:); % First side.
        sideII=boundPoly(sideNum+2,:)-boundPoly(sideNum+1,:); % Second side.
    elseif sideNum==numSides
        sideI=boundPoly(sideNum+1,:)-boundPoly(sideNum,:); % Last side.
        sideII=boundPoly(2,:)-boundPoly(1,:); % 1st side again.
    end
    
    sideLengths(sideNum)=norm(sideI);
    angleThetas(sideNum)=180-acosd(dot(sideI,sideII)/(norm(sideI)*norm(sideII))); % Angle between two sides.
    
end

%% 3. Determine angle between sides and positive x axis.
sideThetas=zeros(numSides,1);
for sideNum=1:numSides
    if sideNum<numSides
        sideI=boundPoly(sideNum+1,:)-boundPoly(sideNum,:); % First side.
        %         sideII=boundPoly(sideNum+2,:)-boundPoly(sideNum+1,:); % Second side.
    elseif sideNum==numSides
        sideI=boundPoly(sideNum+1,:)-boundPoly(sideNum,:); % Last side.
        %         sideII=boundPoly(2,:)-boundPoly(1,:); % 1st side again.
    end
    
    sideThetas(sideNum)=atan2d(sideI(2),sideI(1)); % Angle between side and positive horizontal axis.
    
end

%% 4. Rotate the convex hull (such that the first side of widest angle is on x ais).
if isequal(angleThetas(1)*ones(size(angleThetas)),angleThetas) % If all the interior angles are the same.
    maxAngleTheta=angleThetas(1);
    thetaIdx=1;
else % If there is a maximum.
    [maxAngleTheta,~]=max(angleThetas(angleThetas<180 & angleThetas>0)); % e.g. angleNum=3 means between side 3 & 4. angleNum=5 means between side 5 & 1.
    thetaIdx=find(angleThetas==maxAngleTheta,1);
end
theta=-1*sideThetas(thetaIdx); % Amount to rotate.
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)]; % Rotation matrix.
rotPoly = R * boundPoly';
rotPoly=rotPoly';

% Check if rotated polygon contains vertical lines.
xdiffs=diff(rotPoly(:,1));
if any(xdiffs) % Check if two neighboring x values are the same (vertical line).
    
end

% How to rotate the data back to original.
backTheta=theta*-1;
Rback = [cosd(backTheta) -sind(backTheta); sind(backTheta) cosd(backTheta)]; % Rotation matrix.

if plotRotOn==1
    if plotOrigOn~=1
        Q=figure;
    end
    hold on;
    plot(rotPoly(:,1),rotPoly(:,2),'b');
    axis equal;
end

%% 5. Determine slopes of each side of rotPoly
rotSlopes=zeros(numSides,1);
for i=1:numSides
    
    rotSlopes(i)=(rotPoly(i+1,2)-rotPoly(i,2))/(rotPoly(i+1,1)-rotPoly(i,1));
    
end

%% 6. Isolate vertex of two sides with widest angle.
balloonNum=1;

% Rotated polygon, starting from vertex point of widest angle.
initVert=rotPoly(thetaIdx+1,:); % Vertex of widest angle. Is thetaIdx never the last value? 

%% 7. Create bisecting vector from vertex to arbitrary initial interior point.
initDist=0.3*min(sideLengths); % Hypotenuse (distance)
% NEED TO AVOID INITIALIZING THE LENGTH OF THE BISECTING VECTOR SO LONG THAT THE INITIAL POINT IS OUTSIDE OF THE POLYGON.

% % Rotated polygon. Half the widest angle (from horizontal to angle of next side).
halfTheta=180-maxAngleTheta/2; % Angle (degrees) from positive x axis for bisecting vector.
xDir=cosd(halfTheta); % Negative if > 90 degrees. Is the X component of initial bisecting vector
yDir=sind(halfTheta); % Y component of initial bisecting vector.
% if abs(halfTheta)>90 % Need to invert the triangle to determine point coordinates. (bisecting vector in 2nd quadrant)
%     triTheta=180-halfTheta;
%     xDist=-1*cosd(triTheta)*initDist;
% else % Don't need to invert the triangle to determine point coordinates. (bisecting vector in quadrant 1)
%     triTheta=halfTheta;
%     xDist=cosd(triTheta)*initDist;
% end
% m1=yDir/xDir; % Slope of the bisecting vector
% For m2, need to ensure it's not one of the sides containing the initVert.
% corner1=rotPoly(rotPoly(:,1)==min(rotPoly(:,1)),1:2);
% corner1=corner1(1,:); % If point is doubled up, only take the first one.
% corner2=rotPoly(rotPoly(:,2)==max(rotPoly(:,2)),1:2);
% corner2=corner2(1,:); % If point is doubled up, only take the first one.
% if thetaIdx>1
%     cornerNumToUse=thetaIdx-1; % To get the m2 slope. This is just any side that is not one of the existing sides of the initVert.
% else
%     cornerNumToUse=size(rotPoly,1); % To get the m2 slope. This is just any side that is not one of the existing sides of the initVert.
% end
% if cornerNumToUse>1
%     cornerNum2ToUse=cornerNumToUse-1;
% else
%    cornerNum2ToUse=size(rotPoly,1)-1; 
% end
% meanCorners=mean([corner1; corner2],1);
% m2=(corner2(2)-corner1(2))/(corner2(1)-corner1(1)); % May not align with existing sides.
% b1=initVert(2)-m1*initVert(1); % Y-intercept of bisecting vector.
% b2=corner1(2)-m2*corner1(1); % Y-intercept of side 2.
% bothSlopes=(m1-m2); % 1 - 2. X's on left side.
% bothBs=(b2-b1); % 2 - 1. B's on R side.

% commX=bothBs/bothSlopes; % X value of side intersections.
% commY=m1*(commX)+b1; % Y value of side intersections.
if thetaIdx<size(rotPoly,1)-1 % Get the two sides of the initVert.
    initCloseSideNums=[thetaIdx thetaIdx+1];    
else
    initCloseSideNums=[1 thetaIdx];
end
in=0; iterNum=0;
while in==0 % Checks to see if the initial point is inside or outside of the polygon.
    iterNum=iterNum+1;
%     initPoint=initVert+[xDir yDir]*initDist*norm([commX commY]-initVert);
    initPoint=initVert+[xDir yDir]*initDist;
    [in,on]=inpolygon(initPoint(1),initPoint(2),rotPoly(:,1),rotPoly(:,2));
    if ~in || on % Not inside, or is hitting the boundary line.
        initDist=initDist/2; % Keeps shrinking the scaling factor until
    else
        [~,closestSideNum]=findClosestSide(rotPoly,initPoint(1),initPoint(2));
        if ~ismember(closestSideNum,initCloseSideNums) % Ensures that the initial point is closer to the sides of its vertex than other sides (with sharp triangles).
            in=0;
            initDist=initDist/2;
        end
    end
end
% initPoint=initVert+(meanCorners-initVert)/2;
bisectVect{balloonNum}=initPoint-initVert; % From vertex to interior point.

th=0:pi/50:2*pi; % Creates a circle for plotting.
if plotPrevBals==1 || plotRotOn==1
    [dist1,~]=findClosestSide(rotPoly,initPoint(1),initPoint(2));
    xun{balloonNum}=dist1*cos(th)+initPoint(1); % 1st side circle
    yun{balloonNum}=dist1*sin(th)+initPoint(2);
end

if plotPrevBals==1
    % Plot.
    if plotRotOn==1
        scatter(initPoint(1),initPoint(2),'r*');
        plot(xun{balloonNum},yun{balloonNum},'r*');
    end
    if plotOrigOn==1
        backInitPoint=Rback*initPoint';
        backun1=Rback*[xun{balloonNum}; yun{balloonNum}];
        scatter(backInitPoint(1),backInitPoint(2),'r*');
        plot(backun1(1,:),backun1(2,:),'r*');
    end
end

%% 8. Increment along bisecting vector until the next closest side is hit.
newPoint{balloonNum}=initPoint;
iterNum=0;

% [~,initCloseSideNums(1),initCloseSideNums(2)]=findClosestSide(rotPoly,newPoint{balloonNum}(1),newPoint{balloonNum}(2));

a=0;
scaleFactor=0.001*min(sideLengths);
while a==0
    iterNum=iterNum+1;
    newPoint{balloonNum}=newPoint{balloonNum}+bisectVect{balloonNum}*iterNum*scaleFactor; % Move the interior point along the bisecting vector.
    
    [perpDist{balloonNum},newCloseSideNum{balloonNum}]=findClosestSide(rotPoly,newPoint{balloonNum}(1),newPoint{balloonNum}(2));
    
    if ~ismember(newCloseSideNum{balloonNum},initCloseSideNums) % If a new side is contacted.
        if plotPrevBals==1 || plotRotOn==1 || plotOrigOn==1
            xun{balloonNum}=perpDist{balloonNum}*cos(th)+newPoint{balloonNum}(1); % 1st side circle
            yun{balloonNum}=perpDist{balloonNum}*sin(th)+newPoint{balloonNum}(2);
        end
        if plotPrevBals==1
            if plotRotOn==1
                scatter(newPoint{balloonNum}(1),newPoint{balloonNum}(2),'g*');
                plot(xun{balloonNum},yun{balloonNum},'g*');
            end
            if plotOrigOn==1
                backNewPoint=Rback*newPoint{balloonNum}';
                backun2=Rback*[xun{balloonNum}; yun{balloonNum}];
                scatter(backNewPoint(1),backNewPoint(2),'g*');
                plot(backun2(1,:),backun2(2,:),'g*');
            end
        end
        % Stop processing if new side and either of previous sides are parallel.
        if isequal(rotSlopes(newCloseSideNum{balloonNum}),rotSlopes(initCloseSideNums(1))) ...
                || isequal(rotSlopes(newCloseSideNum{balloonNum}),rotSlopes(initCloseSideNums(2)))
            maxRadOrigin=Rback*newPoint{balloonNum}';
            maxRad=perpDist{balloonNum};
            return;
        end
        break;
    end
    
end

%% 9. Determine which of the first two sides to use (first balloon only).
% Use either the first or second side of the first angle
% Whichever side creates the smallest angle with the new side.

if rotSlopes(newCloseSideNum{balloonNum})>0 || isequal(abs(rotSlopes(newCloseSideNum{balloonNum})),inf)
    % Use second side.
    newCloseSides(1:balloonNum+1)=[initCloseSideNums(2) newCloseSideNum{balloonNum}];
elseif rotSlopes(newCloseSideNum{balloonNum})<0
    % Use first side (slope=0).
    newCloseSides(1:balloonNum+1)=[initCloseSideNums(1) newCloseSideNum{balloonNum}];
elseif isequal(rotSlopes(newCloseSideNum{balloonNum}),0)
    maxRadOrigin=Rback*newPoint{balloonNum}';
    maxRad=perpDist{balloonNum};
    return;
end

%% 10. While loop to iterate through balloons.
% 1. Have 2 new sides to find intersection of. (newCloseSides var)
% 2. Find intersection.
% 3. Compute unit vector between center of current balloon and vertex (in direction of balloon).
% 4. Iterate until a third side is hit.
% 5. Record the 2 newest side numbers.
% 6. Repeat 1-5 UNTIL step 4 results in a smaller distance than previous distance.

% Iterate through balloons.
while a==0
    balloonNum=balloonNum+1;
    newPoint{balloonNum}=newPoint{balloonNum-1};
    
    %% 11. Find where the two newly selected sides intersect (to create bisect vector).
    %     newCloseSides{balloonNum}=sort(newCloseSides{balloonNum});
    [newVert{balloonNum}(1),newVert{balloonNum}(2)]=SideVertFind(rotPoly,newCloseSides(balloonNum-1,1),newCloseSides(balloonNum-1,2));
    if isequal([Inf Inf],abs(newVert{balloonNum})) || any(isnan(newVert{balloonNum}))
        maxRad=perpDist{balloonNum-1};
        maxRadOrigin=Rback*newPoint{balloonNum-1}';
        break;
    end
    
    newBisectVect{balloonNum}=(newPoint{balloonNum}-newVert{balloonNum})/norm(newPoint{balloonNum}-newVert{balloonNum}); % Unit vector.
    
    %% 12. Find next closest side.
    iterNum=0;
    while a==0
        iterNum=iterNum+1;
        newPoint{balloonNum}=newPoint{balloonNum}+newBisectVect{balloonNum}*iterNum*scaleFactor; % Move the interior point along the bisecting vector.
        
        [perpDist{balloonNum},newCloseSideNum{balloonNum}]=findClosestSide(rotPoly,newPoint{balloonNum}(1),newPoint{balloonNum}(2));
        % Termination condition for entire process.
        if abs(perpDist{balloonNum})<abs(perpDist{balloonNum-1})
            maxRad=perpDist{balloonNum-1};
            maxRadOrigin=Rback*newPoint{balloonNum-1}';
            break;
        end
        
        if ~ismember(newCloseSideNum{balloonNum},newCloseSides(balloonNum-1,:)) % If a new side is contacted.
            if plotPrevBals==1 || plotRotOn==1 || plotOrigOn==1
                xun{balloonNum}=perpDist{balloonNum}*cos(th)+newPoint{balloonNum}(1); % 1st side circle
                yun{balloonNum}=perpDist{balloonNum}*sin(th)+newPoint{balloonNum}(2);
            end
            if plotPrevBals==1
                if plotRotOn==1
                    scatter(newPoint{balloonNum}(1),newPoint{balloonNum}(2),'c*');
                    plot(xun{balloonNum},yun{balloonNum},'c*');
                end
                if plotOrigOn==1
                    backNewPoint2=Rback*newPoint{balloonNum}';
                    backun3=Rback*[xun{balloonNum}; yun{balloonNum}];
                    scatter(backNewPoint2(1),backNewPoint2(2),'c*');
                    plot(backun3(1,:),backun3(2,:),'c*');
                end
            end
            % Generate the two new sides for the next balloon. Keep new one, remove old one.
            newCloseSides(balloonNum,1)=newCloseSides(balloonNum-1,2);
            newCloseSides(balloonNum,2)=newCloseSideNum{balloonNum};
            % If two newest sides are parallel, stop all execution here.
            break;
        else % No new side is hit
            % Current sides are parallel.
            if isequal(sideThetas(newCloseSides(balloonNum-1,1)),sideThetas(newCloseSides(balloonNum-1,2))) ...
                    || isequal(round(abs(sideThetas(newCloseSides(balloonNum-1,1)))+abs(sideThetas(newCloseSides(balloonNum-1,2)))),180)
                maxRad=perpDist{balloonNum-1};
                maxRadOrigin=Rback*newPoint{balloonNum-1}';
                return;
            end
        end
        
    end
    
    % Termination condition for entire process.
    if abs(perpDist{balloonNum})<abs(perpDist{balloonNum-1})
        if plotRotOn==1
            scatter(newPoint{balloonNum-1}(1),newPoint{balloonNum-1}(2),'b*');
            plot(xun{balloonNum-1},yun{balloonNum-1},'b*');
        end
        break;
    end
    
end

% disp(maxRad);

%% 11. Plot de-rotated location.
% disp(maxRadOrigin);
if plotOrigOn==1
    scatter(maxRadOrigin(1),maxRadOrigin(2),'k*');
    xunFin=perpDist{balloonNum-1}*cos(th)+maxRadOrigin(1); % 1st side circle
    yunFin=perpDist{balloonNum-1}*sin(th)+maxRadOrigin(2);
    plot(xunFin,yunFin,'k*');
end

if plotOrigOn==1 || plotRotOn==1
    close(Q);
end

end

function [minDist,closestSideNum]=findClosestSide(rotPoly,xP,yP)
% Returns the side number of the closest side to an interior point of convex hull.
x=[xP yP];
for aa=1:size(rotPoly,1)-1 % Iterate through points.
    %     sides(a,1:2)=rotPoly(a+1,:)-rotPoly(a,:);
    a=rotPoly(aa,:);
    b=rotPoly(aa+1,:);
    
    d_ab = norm(a-b);
    d_ax = norm(a-x);
    d_bx = norm(b-x);
    
    if dot(a-b,x-b)*dot(b-a,x-a)>=0
        A = [a,1;b,1;x,1];
        dists(aa) = abs(det(A))/d_ab;
    else
        dists(aa) = min(d_ax, d_bx);
    end
    
end
[minDist,closestSideNum]=min(abs(dists));
end

function [commX,commY]=SideVertFind(boundPoly,row_sideNum1_vert,row_sideNum2_vert)

%% Inputs: 
% boundPoly: The boundary of the polygon (from Matlab boundary function)
% row_sideNum1_vert: The row number in the boundPoly variable of the first side's first vertex.
% row_sideNum2_vert: The row number in the boundPoly variable of the second side's first vertex.
%% Outputs:
% (commX,commY): Common X and Y values of both sides (point of intersection).

sideSlope1=(boundPoly(row_sideNum1_vert+1,2)-boundPoly(row_sideNum1_vert,2))/(boundPoly(row_sideNum1_vert+1,1)-boundPoly(row_sideNum1_vert,1));
sideSlope2=(boundPoly(row_sideNum2_vert+1,2)-boundPoly(row_sideNum2_vert,2))/(boundPoly(row_sideNum2_vert+1,1)-boundPoly(row_sideNum2_vert,1));

b1=boundPoly(row_sideNum1_vert,2)-sideSlope1*boundPoly(row_sideNum1_vert,1); % Y-intercept of side 1.
b2=boundPoly(row_sideNum2_vert,2)-sideSlope2*boundPoly(row_sideNum2_vert,1); % Y-intercept of side 2.
bothSlopes=(sideSlope1-sideSlope2); % 1 - 2. X's on left side.
bothBs=(b2-b1); % 2 - 1. B's on R side.
commX=bothBs/bothSlopes; % X value of side intersections.
commY=sideSlope1*(commX)+b1; % Y value of side intersections.

end