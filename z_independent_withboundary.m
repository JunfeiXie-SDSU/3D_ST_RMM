%=========================================================================
% z-independent smooth-turn random mobility model 
%
%Reference: J. Xie, Y. Wan, B. Wang, S. Fu, K. Lu, J. H. Kim, "A Comprehensive 
%           3-Dimensional Random Mobility Modeling Framework for Airborne Networks",
%           IEEE Access, Vol. 6, pp. 22849-22862.
%
%Date: 03/10/2019
%=========================================================================

%=========================================================================
% z-independent random mobility model demo 

%%%Functionality:
    % Generate a trajectory for fixed-wing aircraft using the z-independent 
    % smooth-turn random mobility model. This model is suitable for
    % applications with small variation in aerial mobility along the z
    % direction, such as civilian and commercial applications.

%%%Parameters to be configured:
    % deltaT: sampling rate (s)
    % t_simulation: total simulation time (s), i.e.,length of the trajectory
    % (x,y,z):  initial position vector (m) of the aircraft
    % alpha, vxy_mu, vxy_sig: parameters used to calculate the speed of the aircraft along the x-y plane according to the equation: v=alpha*v+(1-alpha)*vxy_mu+sqrt((1-alpha^2))*vxy_sig*randn;
    % az_mu: mean vertical acceleration (m/s^2), where the vertical acceleration is assumed to follow a truncated exponential distribution
    % vz_min: minimum vertical speed (m/s)
    % vz_max: maximum vertical speed (m/s)
    % az_min: minimum vertical acceleration (m^2/s)
    % az_max: maximum vertical acceleration (m^2/s)
    % v_xy,v_z: initial speed of the aircraft (m/s) along the x-y plane and z direction, respectively.
    % R_mu: mean value of the turn radius (m), whose inverse absolute value, i.e., 1/|R|, is assumed to follow a truncated exponential distribution
    % R_min: minimum turn radius (m)
    % T_mu: mean value of the travel time interval (s), which is assumed to follow an exponential distribution
    % R_positive_proba: probability of positive turn radius, where positive (negative) turn radius represents right (left) turns  
    % az_positive_proba: probability of positive vertical acceleration
    % phi: heading angle of the aircraft measured anti-clockwise
    % theta: turn angle of the aircraft
    % x_boundary,y_boundary,z_boundary: define the boundary of the flight zone
    % x,y,z: initial positions of the aircraft

%=========================================================================

clc;
clear all;

% Configure parameters
deltaT=0.1;
t_simulation=200;

alpha=0.94;vxy_mu=15.34;vxy_sig=2.35;
az_mu=0.59;vz_min=-0.5;vz_max=0.5;v_xy=10;vz=.12;az_tmin=-1.0;az_tmax=1.0;
R_mu=20;
R_min=10;
T_mu=0.1987;
R_positive_proba=0.5;
az_positive_proba=0.5;
phi=pi;
theta=0;
x_boundary=100;y_boundary=100;z_boundary=50;
x=x_boundary/2;y=y_boundary/2;z=z_boundary/2;

%Find the buffer zone
dx_boundary=(1+sqrt(1./2))*R_min;
dy_boundary=(1+sqrt(1./2))*R_min;
dz_boundary=vz_max^2./(2*az_tmax);
x_buffer_up=x_boundary-dx_boundary;y_buffer_up=y_boundary-dy_boundary;z_buffer_up=z_boundary-dz_boundary;
x_buffer_down=0+dx_boundary;y_buffer_down=0+dy_boundary;z_buffer_down=0+dz_boundary;

%Draw the boundary
figure;
axis([0 x_boundary 0 y_boundary 0 z_boundary]);
hold on;

while t_simulation>0
    % Randomly generate a travel time interval, which follows an  exponential distribution 
    travel_time_interval=random('exponential',T_mu); 
    % Randomly generate the turn radius, which follows a truncated exponential distribution 
    r=(2*(rand<R_positive_proba)-1)*1./randomNumTrunExpDis(1/R_mu,0,1/R_min,1);
    % Calculate the turn rate 
    W=v_xy/r;
    % Calculate the speed of the aircraft along the x-y plane
    v_xy=alpha*v_xy+(1-alpha)*vxy_mu+sqrt((1-alpha^2))*vxy_sig*randn;   
    % Calculate the position of the turning center
    cx=x+r*sin(phi);
    cy=y-r*cos(phi);
    % Begin to move
    while travel_time_interval>0  
        theta=W*deltaT;
        phi=phi-theta; 
        x=cx-r*sin(phi);
        y=cy+r*cos(phi);
        % Randomly generate the vertical acceleration, which follows a truncated exponential distribution
        az=(2*(rand<az_positive_proba)-1)*randomNumTrunExpDis(az_mu,0,abs(az_tmax),1);
        z=z+vz*deltaT+0.5*az*deltaT^2;
        vz=vz+az*deltaT;
        % Set the vertical speed limit of the aircraft
        if(vz>vz_max) 
            vz=vz_max;
        end
        if(vz<vz_min)
            vz=vz_min;
        end
        plot3(x,y,z,'Marker','o','MarkerSize',3,'MarkerEdgeColor','r');
        hold on;
        xlabel('X(m)','FontSize',16);
        ylabel('Y(m)','FontSize',16);
        zlabel('Z(m)','FontSize',16);
        M=getframe;
        travel_time_interval=travel_time_interval-deltaT;
        t_simulation=t_simulation-deltaT;
        
        %When the aircraft enters the buffer zone
        firstTime=1; % Indicate if the aircraft just enters the buffer zone
        left=0; % Indicate if the aircraft turns left or turns right
        while(x<x_buffer_down||y<y_buffer_down||z<z_buffer_down||x>x_buffer_up||y>y_buffer_up||z>z_buffer_up)
            travel_time_interval=0; %When the aircraft is in buffer zone,  set the travel_time_interval to 0
            %When the aircraft is in the buffer zone along x, y, and z directions
            if((x<x_buffer_down||y<y_buffer_down||x>x_buffer_up||y>y_buffer_up)&&(z<z_buffer_down||z>z_buffer_up)) 
                r=R_min; %Set the turn radius to the minimum turn radius.
                %Determine to turn clockwise or anti clockwise
                [cx,cy,left,firstTime]= determineTurningCenter(x,y,phi,r,x_buffer_down,x_buffer_up,y_buffer_down,y_buffer_up,v_xy,cx,cy,left,firstTime);
                %Determine the turn rate
                W=determineW(left, v_xy,r);
                theta=W*deltaT;
                phi=phi-theta; 
                [x,y]=calculatePosition(left, v_xy,cx,cy,r,phi);
                %When the aircraft enters the bottom of the buffer zone
                if(z<z_buffer_down)
                    %Calculate the acceleration in z direction
                    az_dz=(2*dz_boundary./(vz^2));
                    az = (az_tmax-az_dz).*rand(1,1) + az_dz;
                    z=z+vz*deltaT+0.5*az*deltaT^2;
                    vz=vz+az*deltaT;
                    if(vz>vz_max)
                        vz=vz_max;
                    end
                    if(vz<vz_min)
                        vz=vz_min;
                    end
                end
                
                %When the aircraft enters the top of the buffer zone
                if(z>z_buffer_up)
                    %Calculate the acceleration in z direction
                    az_dz=(2*dz_boundary./(vz^2));
                    az = -((az_tmax-az_dz).*rand(1,1) + az_dz); 
                    z=z+vz*deltaT+0.5*az*deltaT^2;
                    vz=vz+az*deltaT;
                    if(vz>vz_max)
                        vz=vz_max;
                    end
                    if(vz<vz_min)
                        vz=vz_min;
                    end
                end
                plot3(x,y,z,'Marker','o','MarkerSize',3,'MarkerEdgeColor','r');
                hold on;
                xlabel('X(m)','FontSize',16);
                ylabel('Y(m)','FontSize',16);
                zlabel('Z(m)','FontSize',16);
                M=getframe;
                t_simulation=t_simulation-deltaT;
            end 
            
            %When the aircraft is in the buffer zone along x and y directions
            if((x<x_buffer_down||y<y_buffer_down||x>x_buffer_up||y>y_buffer_up)&&(z>=z_buffer_down&&z<=z_buffer_up))  
                r=R_min;
                [cx,cy,left,firstTime]= determineTurningCenter(x,y,phi,r,x_buffer_down,x_buffer_up,y_buffer_down,y_buffer_up,v_xy,cx,cy,left,firstTime);
                W=determineW(left, v_xy,r);
                theta=W*deltaT;
                phi=phi-theta; 
                [x,y]=calculatePosition(left, v_xy,cx,cy,r,phi);
                az=(2*(rand<az_positive_proba)-1)*randomNumTrunExpDis(az_mu,0,abs(az_tmax),1); 
                z=z+vz*deltaT+0.5*az*deltaT^2;
                vz=vz+az*deltaT;
                if(vz>vz_max)
                   vz=vz_max;
                end
                if(vz<vz_min)
                   vz=vz_min;
                end
                plot3(x,y,z,'Marker','o','MarkerSize',3,'MarkerEdgeColor','r');
                hold on;
                xlabel('X(m)','FontSize',16);
                ylabel('Y(m)','FontSize',16);
                zlabel('Z(m)','FontSize',16);
                M=getframe;
                t_simulation=t_simulation-deltaT;
            end 
            
            %When the aircraft is in the buffer zone along the z direction
            if((x>=x_buffer_down&&y>=y_buffer_down&&x<=x_buffer_up&&y<=y_buffer_up)&&(z<z_buffer_down||z>z_buffer_up)) 
                theta=W*deltaT;
                phi=phi-theta; 
                x=cx-r*sin(phi);
                y=cy+r*cos(phi);
              
                if(z<z_buffer_down)
                    az = (az_tmax-(2*dz_boundary./(vz^2))).*rand(1,1) + (2*dz_boundary./(vz^2));
                    z=z+vz*deltaT+0.5*az*deltaT^2;
                    vz=vz+az*deltaT;
                    if(vz>vz_max)
                        vz=vz_max;
                    end
                    if(vz<vz_min)
                        vz=vz_min;
                    end
                end
                
                if(z>z_buffer_up)
                    az = -((az_tmax-(2*dz_boundary./(vz^2))).*rand(1,1) + (2*dz_boundary./(vz^2)));
                    z=z+vz*deltaT+0.5*az*deltaT^2;
                    vz=vz+az*deltaT;
                    if(vz>vz_max)
                        vz=vz_max;
                    end
                    if(vz<vz_min)
                        vz=vz_min;
                    end
                end
                plot3(x,y,z,'Marker','o','MarkerSize',3,'MarkerEdgeColor','r');
                hold on;
                xlabel('X(m)','FontSize',16);
                ylabel('Y(m)','FontSize',16);
                zlabel('Z(m)','FontSize',16);
                M=getframe;
                travel_time_interval=travel_time_interval-deltaT;
                t_simulation=t_simulation-deltaT;
            end      
        end
    end
end


%Determine if the aircraft is in the buffer zone along x and y directions
function result=isXYOutBoundary(x,y,x_buffer_down,x_buffer_up,y_buffer_down,y_buffer_up)
    if(x<x_buffer_down||y<y_buffer_down||x>x_buffer_up||y>y_buffer_up)
        result=1;
    else
        result=0;
    end
end


%Determine the turning center of the aircraft
function [cx,cy,left,firstTime]= determineTurningCenter(x,y,phi,r,x_buffer_down,x_buffer_up,y_buffer_down,y_buffer_up,v_xy,icx,icy,ileft,firstTime)
    if(~firstTime)
       cx=icx;
       cy=icy;
       left=ileft;
    else
        if(v_xy>0)
            cx_right=x+r*sin(phi);
            cy_right=y-r*cos(phi);
            cx_left=x-r*sin(phi);
            cy_left=y+r*cos(phi);
        else
            cx_left=x+r*sin(phi);
            cy_left=y-r*cos(phi);
            cx_right=x-r*sin(phi);
            cy_right=y+r*cos(phi);
        end

        cx=cx_left;
        cy=cy_left;
        left=1;
        if(isXYOutBoundary(cx_right,cy_right,x_buffer_down,x_buffer_up,y_buffer_down,y_buffer_up)&&~isXYOutBoundary(cx_left,cy_left,x_buffer_down,x_buffer_up,y_buffer_down,y_buffer_up))
            cx=cx_left;
            cy=cy_left;
            left=1;
        end

        if(~isXYOutBoundary(cx_right,cy_right,x_buffer_down,x_buffer_up,y_buffer_down,y_buffer_up)&&isXYOutBoundary(cx_left,cy_left,x_buffer_down,x_buffer_up,y_buffer_down,y_buffer_up))
            cx=cx_right;
            cy=cy_right;
            left=0;
        end

        if(isXYOutBoundary(cx_right,cy_right,x_buffer_down,x_buffer_up,y_buffer_down,y_buffer_up)&&isXYOutBoundary(cx_left,cy_left,x_buffer_down,x_buffer_up,y_buffer_down,y_buffer_up))
            if(cx_right<=x_buffer_up&&cx_right>=x_buffer_down&&cy_right>y_buffer_up&&cx_left<x_buffer_down&&cy_left<=y_buffer_up&&cy_left>=y_buffer_down)

                if((cy_right-y_buffer_up)>(x_buffer_down-cx_left))
                    cx=cx_left;
                    cy=cy_left;
                    left=1;
                else
                    cx=cx_right;
                    cy=cy_right;
                    left=0;
                end
            end

            if(cx_left<=x_buffer_up&&cx_left>=x_buffer_down&&cy_left>y_buffer_up&&cx_right>x_buffer_up&&cy_right<=y_buffer_up&&cy_right>=y_buffer_down)

                if((cx_right-x_buffer_up)>(cy_left-y_buffer_up))
                    cx=cx_left;
                    cy=cy_left;
                    left=1;
                else
                    cx=cx_right;
                    cy=cy_right;
                    left=0;
                end
            end

            if(cx_left>x_buffer_up&&cy_left>=y_buffer_down&&cy_left<=y_buffer_up&&    cx_right<=x_buffer_up&&cx_right>=x_buffer_down&&cy_right<y_buffer_down)        
                if((y_buffer_down-cy_right)>(cx_left-x_buffer_up))
                    cx=cx_left;
                    cy=cy_left;
                    left=1;
                else
                    cx=cx_right;
                    cy=cy_right;
                    left=0;
                end
            end
            if(cx_right<x_buffer_down&&cy_right>=y_buffer_down&&cy_right<=y_buffer_up&&    cx_left<=x_buffer_up&&cx_left>=x_buffer_down&&cy_left<y_buffer_down)        
                if((x_buffer_down-cx_right)>(y_buffer_down-cy_left))
                    cx=cx_left;
                    cy=cy_left;
                    left=1;
                else
                    cx=cx_right;
                    cy=cy_right;
                    left=0;
                end
            end 
        end
        firstTime=0;
    end
end


%Calculate the position of the aircraft
function [x,y]=calculatePosition(left, v_xy,cx,cy,r,phi)  
    if(v_xy>0)
        if(left)
            x=cx+r*sin(phi);
            y=cy-r*cos(phi);
        else
            x=cx-r*sin(phi);
            y=cy+r*cos(phi);
        end
    else
        if(~left)
            x=cx+r*sin(phi);
            y=cy-r*cos(phi);
        else
            x=cx-r*sin(phi);
            y=cy+r*cos(phi);
        end
    end
end


%Determine the turn rate of the aircraft
function W=determineW(left, v_xy,r)  
    if(v_xy>0)
        if(left)
            W=-v_xy/r;
        else
            W=v_xy/r;
        end
    else
        if(~left)
            W=-v_xy/r;
        else
            W=v_xy/r;
        end
    end
end


% Generate a random number from a truncated exponential distribution
function result=randomNumTrunExpDis(mu,lower_limit,upper_limit,timeLength)
    time=0*mu:mu/100:2*mu;
    len=length(time);
    pd=makedist('Exponential','mu',mu);
    pd=truncate(pd,lower_limit,upper_limit);
    value=[];

    for i=1:len
        temp=cdf(pd,time(i));
        if(temp>0.999)
            temp=1;
        end
        value=[value,temp];
    end

    [value, mask] = unique(value);
    time= time(mask);
    randomValues=rand(1,timeLength);
    while(randomValues>max(value))
        randomValues=rand(1,timeLength);
    end
    result=interp1(value,time,randomValues);
end