%=========================================================================
% z-dependent smooth-turn random mobility model 
%
%Reference: J. Xie, Y. Wan, B. Wang, S. Fu, K. Lu, J. H. Kim, "A Comprehensive 
%           3-Dimensional Random Mobility Modeling Framework for Airborne Networks",
%           IEEE Access, Vol. 6, pp. 22849-22862.
%
%Date: 03/10/2019
%=========================================================================

%=========================================================================
% z-dependent random mobility model demo 

%%%Functionality:
    % Generate a trajectory for fixed-wing aircraft using the z-dependent 
    % smooth-turn random mobility model. This model is suitable for
    % applications involving climbing or descending turns, such as military
    % applications and air shows. 
    
%%%Parameters to be configured:
    % deltaT: sampling rate (s)
    % t_simulation: total simulation time (s), i.e.,length of the trajectory
    % alpha, v_mu, v_sig: parameters used to calculate the velocity of the aircraft according to the equation: v=alpha*v+(1-alpha)*v_mu+sqrt((1-alpha^2))*v_sig*randn;
    % v: initial velocity vector of the aircraft (m/s)
    % R_mu: mean of the turn radius (m), whose inverse absolute value, i.e., 1/|R|, is assumed to follow a truncated exponential distribution
    % R_min: minimum turn radius (m)
    % T_mu: mean of the travel time interval (s), which is assumed to follow an exponential distribution
    % R_positive_proba: probability of the positive radius
    % x_boundary, y_boundary, z_boundary: define the boundary of the flight zone
    % x, y, z: initial positions of the aircraft
%=========================================================================

clc;
clear all;

% Configure parameters
deltaT=0.1;
t_simulation=200;
alpha=0.94;v_mu=12.06;v_sig=3.04;
v=[10;10;1];
R_mu=100;
R_min=50;
T_mu=0.6349;
R_positive_proba=0.5;

% Set the boundary
x_boundary=1000;y_boundary=1000;z_boundary=500;
x=1/2*x_boundary;
y=1/2*y_boundary;
z=1/2*z_boundary;
p=[x;y;z];

%Find the buffer zone
dx_boundary=2*R_min;
dy_boundary=2*R_min;
dz_boundary=2*R_min;

x_buffer_up=x_boundary-dx_boundary;y_buffer_up=y_boundary-dy_boundary;z_buffer_up=z_boundary-dz_boundary;
x_buffer_down=0+dx_boundary;y_buffer_down=0+dy_boundary;z_buffer_down=0+dz_boundary;

%Draw the boundary
figure;
axis([0 x_boundary 0 y_boundary 0 z_boundary]);
hold on;

while (t_simulation>0)
    % Randomly generate a travel time interval, which follows an  exponential distribution 
    travel_time_interval=random('exponential',T_mu);
    % Randomly generate the turn radius, which follows a truncated exponential distribution 
    r=(2*(rand<R_positive_proba)-1)*1./randomNumTrunExpDis(1/R_mu,0,1/R_min,1);
    % When the aircraft enters the buffer zone, set the turn radius to minimum turn radius. 
    if(x<x_buffer_down||y<x_buffer_down||z<x_buffer_down||x_buffer_up||y>x_buffer_up||z>x_buffer_up)
     r=R_min ;
    end
    x=p(1);
    y=p(2);
    z=p(3);
    v=alpha*v+(1-alpha)*v_mu+sqrt((1-alpha^2))*v_sig*randn;
    vx=v(1);
    vy=v(2);
    vz=v(3);
    theta=cos(asin(abs(vz)/sqrt(vx^2+vy^2+vz^2)));
    beta=r*theta;
    [cx, cy,cz]= calculateTurningCenter(p,beta,r,v);
    pc=[cx-p(1);cy-p(2);cz-p(3)];

    % Calculate the acceleration
    at=0;
    an=sum(v.^2)/r;
    aa=sqrt(an^2+at^2);
    a=an*pc./norm(pc)+at*v./norm(v);
    ax=a(1);
    ay=a(2);
    az=a(3);
    
    % Begin to move  
    state=[x;vx;ax;y;vy;ay;z;vz;az];
    w=norm(cross(v,a)./(norm(v)^2));
    A=[1,sin(w*deltaT)/w,(1-cos(w*deltaT))/(w^2);
        0,cos(w*deltaT),sin(w*deltaT)/w;
        0,-w*sin(w*deltaT),cos(w*deltaT)];
    I=[0 0 0;0 0 0;0 0 0];
    

    %When the aircraft is in the buffer zone
    while(state(1)<x_buffer_down||state(4)<y_buffer_down||state(7)<z_buffer_down||state(1)>x_buffer_up||state(4)>y_buffer_up||state(7)>z_buffer_up)
        travel_time_interval=0;
        vi=norm([state(2),state(5),state(8)]);
        state=[A I I;I A I;I I A]*state;
        t_simulation=t_simulation-deltaT;
        plot3(state(1),state(4),state(7),'Marker','o','MarkerSize',3,'MarkerEdgeColor','r');
        hold on;
        xlabel('X(m)','FontSize',16);
        ylabel('Y(m)','FontSize',16);
        zlabel('Z(m)','FontSize',16);
        M=getframe;
        t_simulation=t_simulation-deltaT;
    end
    
   %When the aircraft is out of the buffer zone
    while travel_time_interval>0
        vi=norm([state(2),state(5),state(8)]);
        state=[A I I;I A I;I I A]*state;
        travel_time_interval=travel_time_interval-deltaT;
        t_simulation=t_simulation-deltaT;
        plot3(state(1),state(4),state(7),'Marker','o','MarkerSize',3,'MarkerEdgeColor','r');
        hold on;
        xlabel('X(m)','FontSize',16);
        ylabel('Y(m)','FontSize',16);
        zlabel('Z(m)','FontSize',16);
        M=getframe;
        %When the aircraft enters the buffer zone
        if(state(1)<x_buffer_down||state(4)<y_buffer_down||state(7)<z_buffer_down||state(1)>x_buffer_up||state(4)>y_buffer_up||state(7)>z_buffer_up)
            break;
        end
    end
    p=[state(1),state(4),state(7)];
    v=[state(2);state(5);state(8)];
    a=[state(3);state(6);state(9)];

end

% Calculate the position of the turning center 
function [cx, cy, cz]= calculateTurningCenter(p,beta,r,v)

    cz=p(3)-beta+2*beta*rand;
    syms A B;
    pc=[A-p(1);B-p(2);cz-p(3)];
    rpc=sum(pc.^2)-r^2;
    vpc=sum(v.*pc);
    [sol_a,sol_b]=solve(rpc,vpc,A,B);
    sol_a=double(sol_a);
    sol_b=double(sol_b);
    % If there are two solutions, randomly select one
    if numel(sol_a)==2
        number=randi(2,1);
        cx=sol_a(number);
        cy=sol_b(number);            
    else
        cx=sol_a;
        cy=sol_b;
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

