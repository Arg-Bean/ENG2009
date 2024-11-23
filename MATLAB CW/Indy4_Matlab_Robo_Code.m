%% clear all variables clear command window
bdclose all %close all simulink block
close all %close all plots
clear all %clear all variables
clc %clear comman window

%% 3-D Path Tracing With Inverse Kinematics
%% Introduction
% This example shows how to create a simple 3D robot

%% Construct The Robot
% Create a |rigidBodyTree| object and rigid bodies with their
% associated joints. Specify the geometric properties of each rigid body
% and add it to the robot.

%Create Robot
robot = rigidBodyTree('DataFormat','column','MaxNumBodies',7);

%Specify lengths of arms
L1 = 0.2;
L2 = 0.2;
L3 = 0.2;
L4 = 0.2;
L5 = 0.2;
L6 = 0.2;

%Create first body (with a revolute joint) L1
%
%
body = rigidBody('link1');
%Define body
body.Mass=1; %1 kg (default)
body.CenterOfMass=[0 0 0]; % [0 0 0] m (default) | [x y z] vector
body.Inertia=[1 1 1 0 0 0]; % [1 1 1 0 0 0] kg.m2 (default) | [Ixx Iyy Izz Iyz Ixz Ixy] vector

%Defining the first joint (Revolute, y axis)
joint = rigidBodyJoint('joint1', 'revolute');
setFixedTransform(joint,trvec2tform([0 0 0]));
joint.JointAxis = [0 1 0];
body.Joint = joint;
addBody(robot, body, 'base');



%Create second body (with prismatic joint) L2
%
%
body = rigidBody('link2');
%body definition
body.Mass=1; %1 kg (default)
body.CenterOfMass=[0 0 0]; % [0 0 0] m (default) | [x y z] vector
body.Inertia=[1 1 1 0 0 0]; % [1 1 1 0 0 0] kg.m2 (default) | [Ixx Iyy Izz Iyz Ixz Ixy] vector


%Defining the second joint (Prismatic, y axis)
joint = rigidBodyJoint('joint2','prismatic');
%joint.PositionLimits = [0 0.2]
setFixedTransform(joint, trvec2tform([0,L1,0]));
joint.JointAxis = [0 1 0];
body.Joint = joint;
addBody(robot, body, 'link1');



%Create third body (with prismatic joint) L3
%
%
body = rigidBody('link3');
%body definition
body.Mass=1; %1 kg (default)
body.CenterOfMass=[0 0 0]; % [0 0 0] m (default) | [x y z] vector
body.Inertia=[1 1 1 0 0 0]; % [1 1 1 0 0 0] kg.m2 (default) | [Ixx Iyy Izz Iyz Ixz Ixy] vector


%Defining the third joint (Prismatic, z axis)
joint = rigidBodyJoint('joint3','prismatic');
%joint.PositionLimits = [0 0.2]
setFixedTransform(joint, trvec2tform([0,0,L2]));
joint.JointAxis = [0 0 1];
body.Joint = joint;
addBody(robot, body, 'link2');



%Create fourth body (with prismatic joint) L4
%
%
body = rigidBody('link4');
%body definition
body.Mass=1; %1 kg (default)
body.CenterOfMass=[0 0 0]; % [0 0 0] m (default) | [x y z] vector
body.Inertia=[1 1 1 0 0 0]; % [1 1 1 0 0 0] kg.m2 (default) | [Ixx Iyy Izz Iyz Ixz Ixy] vector


%Defining the fourth joint (Prismatic, x axis)
joint = rigidBodyJoint('joint4','prismatic');
%joint.PositionLimits = [0 0.2]
setFixedTransform(joint, trvec2tform([L3,0,0]));
joint.JointAxis = [1 0 0];
body.Joint = joint;
addBody(robot, body, 'link3');



%Create fifth body (with revolute joint) L5
%
%
body = rigidBody('link5');
%body definition
body.Mass=1; %1 kg (default)
body.CenterOfMass=[0 0 0]; % [0 0 0] m (default) | [x y z] vector
body.Inertia=[1 1 1 0 0 0]; % [1 1 1 0 0 0] kg.m2 (default) | [Ixx Iyy Izz Iyz Ixz Ixy] vector


%Defining the fifth joint (Revolute, z axis)
joint = rigidBodyJoint('joint5','revolute');
setFixedTransform(joint, trvec2tform([0,0,L4]));
joint.JointAxis = [0 0 1];
body.Joint = joint;
addBody(robot, body, 'link4');



%Create sixth body (with revolute joint) L6
%
%
body = rigidBody('link6');
%body definition
body.Mass=1; %1 kg (default)
body.CenterOfMass=[0 0 0]; % [0 0 0] m (default) | [x y z] vector
body.Inertia=[1 1 1 0 0 0]; % [1 1 1 0 0 0] kg.m2 (default) | [Ixx Iyy Izz Iyz Ixz Ixy] vector


%Defining the sixth joint (Revolute, y axis)
joint = rigidBodyJoint('joint6','revolute');
setFixedTransform(joint, trvec2tform([0,L5,0]));
joint.JointAxis = [0 1 0];
body.Joint = joint;
addBody(robot, body, 'link5');


%Create seventh body (with fixed joint) L7
%
%
body = rigidBody('tool');
%body definition
body.Mass=1; %1 kg (default)
body.CenterOfMass=[0 0 0]; % [0 0 0] m (default) | [x y z] vector
body.Inertia=[1 1 1 0 0 0]; % [1 1 1 0 0 0] kg.m2 (default) | [Ixx Iyy Izz Iyz Ixz Ixy] vector


%Defining the seventh joint (fixed)
joint = rigidBodyJoint('fix1','fixed');
setFixedTransform(joint, trvec2tform([L6,0,0]));
body.Joint = joint;
addBody(robot, body, 'link6');
%%
% Show details of the robot to validate the input properties. 
% The robot should have two non-fixed joints for the rigid bodies and a fixed body for the end-effector.
% default values:
% robot.Gravity=[0 0 -9.81] %[0 0 0] (default). i.e. no gravity. To change use: robot.Gravity=[0 0 -9.81]

robot %rigidBodyTree with properties. 
showdetails(robot)

figure(1);
show(robot)

figure(Name = ['Interactive GUI'])
gui = interactiveRigidBodyTree(robot,MarkerScaleFactor=0.5);

%figure(1); 
%show(robot) % show the robot config with zero initial conditions


%% open model
simmdl='WS_3DRobot_example_V2' %simulink model file name
open(simmdl) %open simulink model

%% PID control definition
Kp1=100 %8*3
Kd1=10 %0.1*1
Kfilt1=150

Kp2=100 %8*3
Kd2=10 %0.1*1
Kfilt2=150

Kp3=100 %8*3
Kd3=10 %0.1*1
Kfilt3=150

Kp4=100 %8*3
Kd4=10 %0.1*1
Kfilt4=150

Kp5=100 %8*3
Kd5=10 %0.1*1
Kfilt5=150

Kp6=100 %8*3
Kd6=10 %0.1*1
Kfilt6=150

%% Joint torques definition
% T1amp=1; %Joint 1 torque amplitude
% T2amp=1; %joint 2 torque amplitude
% T3amp=1; %joint 3 torque amplitude
% T4amp=1; %joint 4 torque amplitude
% T5amp=1; %joint 5 torque amplitude
% T6amp=1; %joint 6 torque amplitude

%% Joint angles definition: to try for examples: [45 0], [0 45], [90 -90], [-90 90]
% jointangles=[45 0 0] %deg
% angle1ref=jointangles(1) %Desired joint 1 angle (deg)
% angle2ref=jointangles(2) %Desired joint 2 angle (deg)
% angle3ref=jointangles(3) %Desired joint 3 angle (deg)
%angle4ref=jointangles(4) %Desired joint 4 angle (deg)
%angle5ref=jointangles(5) %Desired joint 5 angle (deg)
%angle6ref=jointangles(6) %Desired joint 6 angle (deg)

%% end effector position definition
initPos=[0.5 -0.5 0.1]; %initial position
initialConfig = homeConfiguration(robot); %inital joint angle

% eePos_des_sim=[0 0.6 0]; %desired position. Default [0 0.6 0];
% %targetPosition =  trvec2tform(eePos_des_sim) %desired position in homogeneous transformation

%% Traj generation
eePos_des_traj=[0.0 0.5 -0.5 0.1;
    10.0 1.0 0.5 0.5;
    20.0 0.5 0.5 0.1;
    30.0 0.5 -0.5 0.1]
time_traj=eePos_des_traj(:,1)
eePos_des_traj=eePos_des_traj(:,2:4)

%% simulate model
tsim=30 %simulation time
t_sample=0.01 %0.01 %simulation time step (sample time)

sim(simmdl) %run the simulation of the model

%% save all variables in the workspace
save('simout')


%% simple plot of inputs and outputs
% figure(1)
% plot(t_sim,input_sim(:,1),'b', t_sim,input_sim(:,2),'--r')
% grid on; hold on;
% title('Input torques vs time'); ylabel('torque (N)'); xlabel('time (sec)')
% legend('torque 1',' torque 2')
% 
% figure(2)
% subplot(2,1,1); plot(t_sim,qs(:,1)*180/pi,'b', t_sim,qs_des(:,1)*180/pi,'--r'); grid on; hold on;
% legend('angle 1 (deg)',' angle 1 desired (deg)');
% ylabel('Angles (deg)'); xlabel('time (sec)');
% title('Joint angles vs time');
% subplot(2,1,2); plot(t_sim,qs(:,2)*180/pi,'b', t_sim,qs_des(:,2)*180/pi,'--r'); grid on; hold on;
% legend('angle 2 (deg)',' angle 2 desired(deg)')
% ylabel('Angles (deg)'); xlabel('time (sec)')
% 
% figure(3)
% subplot(2,1,1); plot(t_sim,eePos_sim(:,1),'b', t_sim,eePos_des(:,1),'--r'); grid on; hold on;
% title('End effector x position vs time'); ylabel('position (m)'); xlabel('time (sec)')
% legend('x','x desired')
% ylabel('x position'); xlabel('time (sec)')
% 
% subplot(2,1,2); plot(t_sim,eePos_sim(:,2),'b', t_sim,eePos_des(:,2),'--r'); grid on; hold on;
% title('End effector y position vs time'); ylabel('position (m)'); xlabel('time (sec)')
% legend('y','y desired')
% ylabel('y position'); xlabel('time (sec)')

% subplot(3,1,3); plot(t_sim,eePos_sim(:,3),'b', t_sim,eePos_des(:,3),'--r'); grid on; hold on;
% title('End effector z position vs time'); ylabel('position (m)'); xlabel('time (sec)')
% legend('z','z desired')
% ylabel('z position'); xlabel('time (sec)')



%% Animate The Solution
% Plot the robot for each frame of the solution using that specific robot configuration.
anim = 0; % control to show animation %0=no anim, 1=anim
gifon = 1; % control whether use gif as output %0=no GIF, 1=create GIF file

%% -------------------------------------------------------------------------------------------------
%%% Animation commands - DO NOT CHNAGE THE COMMAND LINES BELOW, UNLESS YOU KNOW WHAT YOU ARE DOING
%%% -------------------------------------------------------------------------------------------------
% Plot the robot for each frame of the solution using that specific robot
% configuration. Also, plot the desired trajectory.

if exist('anim','var') && anim == true
    %%%
    % Show the robot in the first configuration of the trajectory.
    % Adjust the plot to show the 3-D shape is drawn on.
    % Plot the desired vs actual trajectories.
    figure
    show(robot,qs(1,:)');
    %view(2) %2D view
    view(3) %3D view
    ax = gca;
    ax.Projection = 'orthographic';
    hold on
    %plot(points(:,1),points(:,2),'k')
    axis([-1 1 -1 1])

    plot3(eePos_des_traj(:,1),eePos_des_traj(:,2),eePos_des_traj(:,3),'r--') %desired path
    plot3(eePos_sim(:,1),eePos_sim(:,2),eePos_sim(:,3),'b') %actual path

    %%%
    % Set up an object to display the robot trajectory at a fixed rate of 15 frames per second.
    % Show the robot in each configuration from the inverse kinematic solver.
    % Watch as the arm traces the circular trajectory shown.
    count = length(tsim);
    framesPerSecond = 100; %15;
    r = rateControl(framesPerSecond);

    for i = 1:count
        show(robot,qs(i,:)','PreservePlot',false);
        drawnow
        waitfor(r);
    end
end

%%% GIF Plot Command
% Create gif based on the robot animation

if exist('gifon','var') && gifon == true
    figure
    show(robot,qs(1,:)');
    %view(2) %2D view
    view(3) %3D view
    ax = gca;
    ax.Projection = 'orthographic';
    hold on
    %plot(points(:,1),points(:,2),'k')
    axis([-1 1 -1 1])

    plot3(eePos_des_traj(:,1),eePos_des_traj(:,2),eePos_des_traj(:,3),'r--') %desired path
    plot3(eePos_sim(:,1),eePos_sim(:,2),eePos_sim(:,3),'b') %actual path

    count = length(tsim);
    framesPerSecond = 100; %15;
    r = rateControl(framesPerSecond);

    % change the number 50 to change the file duration (higher means longer
    % GIF video).
    for i = round(linspace(1,count,50))
        show(robot,qs(i,:)','PreservePlot',false);
        drawnow
        waitfor(r);
        exportgraphics(gca,'gifout.gif','Append',true);
    end
end
