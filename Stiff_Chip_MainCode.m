%% Oct. 6, 2026; Author: Yang Liu; this model is designed to simulate the cell deformation and final position under magnetic buoyancy force;
%% Only one unit is simulated. the channel is aligned with the center of the magnet; tapered geometry is applied: gradient 31-7um from bottom to the top of the channel;
%%

%% Initialization of model
clear all;
close all;
format long;
global ka T u_0 M_d eta_f F_r beta M_sf Dia_cell Vol_cell E_cell w_c t_c l_c H_chan W_chan nu XE v_ave Par_magnets d_pf ln_sigma_pf d_pb ln_sigma_pb M_db M_sb M_sc F_head_y l_out


%% Basic constants
ka = 1.38e-23;       % Boltzmann constant, unit m^2*kg*s^(-2)*K^(-1)
T = 273+25;         % Room temperature, where experiments were conducted, unit K
u_0 = pi*4e-7; % Permeability of free space, unit Henry/m

%% Geometries of microchannel (rectangular microchannel)
t_c =10e-6/2; % Half thickness (Z-direction) of channel, unit m, to be measured, by spin curve it should be 150um.
w_c = 1150e-6/2; % Half width (Y-direction) of channel, unit m : % stage two channel 1200um
l_c = 31e-6/2; % Half length (X-direction) of channel, unit m, entrance
l_out = 7e-6/2; % Half length (X-direction) of channel, unit m, exit
s_c = t_c*w_c*4; % Crossectional area (Y-Z plane) of the channel, unit m^2
H_chan = 2*t_c;
W_chan = 2*l_c;

%% Parameters of cells or microbeads
Dia_cell=30*1e-6; % diameter of cell
Vol_cell=(4/3)*pi*(Dia_cell/2)^3;
Phi_cell=0; %volume fraction of on-magnetic particle
Par_beads{1}=[Dia_cell ;Vol_cell;Phi_cell;];
E_cell = 0.02; % stiffness of cell, unit: pa
nu=0.49; %~ 0.49 for nearly incompressible cytoplasm.

%% Parameters of ferrofluids (magnetic nanoparticle within and its other fluidic properties)
M_d = 370e+3; % Ferrofluid nanoparticle bulk volumetric magnetization, unit A/m. For bulk magnetization in (u_0*Mb_f) representation: magnetite (Fe3O4) is 5600 Gauss or 0.56 T (Rosensweig book); maghemite (y-Fe2O3) is 4650 Gauss or 0.465 T (Doyle paper). Particle maghemite bulk magnetization was reported to be lower than bulk material, 310-330kA/m, by Fonnum et al (JMMM 293 41 2005) and Johansson et al (JMMM 173 5 1997)
d_pf = 11.16e-9; % Ferrofluid particle mean diameter (log-normal), unit m
ln_sigma_pf = 0.4426; % ln sigma of ferrofluid nanoparticle size distribution (log-normal), unitless
ro = 1060.6; % Ferrofluid density (considering both magnetic cores and surfactants), kg/m^3
Phi_f = 0.3; % Ferrofluid nanoparticle volume fraction (only considering magnetic cores), unitless
eta_f = 1.108e-3*1/(1-127.3*Phi_f+4505*Phi_f^2); % Ferrofluid viscoity (considering magnetic solids alone), unit Pa.s. Fitted eta_f [mPa.s] = 1/(1-2.464*Phi_f+0.8275*Phi_f^2) [%], Rosensweig book. only valid within 0.028% and 0.29%.
%eta_f=0.99e-3;
M_sf = Phi_f*M_d; % Ferrofluid saturation volumetric magnetization, unit A/m


%% Parameter of magnetic nanoparticles in microbeads (Dynabead MyOne)
M_db = 370e+3; % Magnetic nanoparticle bulk magnetization in Dynabead MyOne, unit A/m, assuming to be maghemite with (u_0*Mb_p) of 4650 Gauss or 0.465 T
d_pb = 11.16e-9; % mean magnetic nanoparticle diameter in microbeads
ln_sigma_pb = 0.44; % ln sigma of magnetic nanoparticle diameter in microbeads (log-normal), unitless
Phi_p = 0.3; % volume fraction of magnetic content in microbead volume, unitless
M_sb = Phi_p*M_db; % Magnetic bead saturation volumetric magnetization, unit A/m

%% Parameter of flow pressure due to heigh difference
g = 9.81; % gravity
h = 6e-3; % liquid height difference between inlets/outlets
Delta_p = ro * g * h;
A_ch   = H_chan * W_chan;
A_cell = pi * (Dia_cell/2)^2;
psi= min(1, A_cell / max(A_ch, eps));
F_head_y = Delta_p * A_ch * psi;

%% Parameters of flows or particles
Q=0;%flow rate, unit: ul/min
F_r = Q*1e-9/60; % Total flow rate, m^3/s
v_ave = F_r/s_c; % Average velocity of flow, m/s

%% Pre-processing for ODE solver
Intime = 1000; % Estimated integral time or particle travelling time in channel, s
Intnum = 10000; % Number of intervals within integral, previously set to 10,000
tspan = linspace(0,Intime,Intnum);
XE = w_c; % Integral ending position (at the end of length of the channel w_c), unit m
options = odeset('InitialStep',0.001,'RelTol', 1e-11,'Events',@events_end_of_channel);

%% Parameters of magnets [Remanent magnetization in x y z directions, unit T; dimensions of magent in x y z directions, unit m; center of permament magnet in microchannal coodinates, unit m; rotational angles of permanent magnet wrt x, y, z axis]
% Remanent magnetization of permanent magnets, unit T. Magnetization is
% fixed to be in y direction of magnet coordinates.
Magnet_Ms = 1.54; %1.1574; % magnetization of permanent magnet, unit T
t_m = 6.35e-3/2;  % half thickness (z-direction) of permanent magnet, unit m
w_m = 6.35e-3/2;  % half width (y-direction) of permanent magnet, unit m
l_m = 38.1e-3/2;  % half length (x-direction) of permanet magnet, unit m
delta_z = 100e-6; % distance between channel wall (close to magnet) and magnet surface (close to channel), unit m
delta_Y = 0e-6; % distance between channel wall (close to magnet) and magnet surface (close to channel), unit m


% magnetization is fixed in y-direction of permanent magnet body coooridnates
Magnet_1_Ms = [0 Magnet_Ms 0];

% Dimenions of permanent magnets, unit m
Magnet_1_Dimension = [l_m*2 w_m*2 t_m*2];


% Center of permanent magnets, relative to center of the microchannel,
% which is always at [0 0 0], unit m
Magnet_1_Center = [0 -w_m-w_c-delta_Y -t_m-t_c-delta_z];

Magnet_1_Angle = [0 0 0];    % repelling type #2


% Creat a structure to store all above infomration for permanent magents
Par_magnets{1} = [Magnet_1_Ms; Magnet_1_Dimension; Magnet_1_Center; Magnet_1_Angle];
%%magnetic field


%% Plot the positions of microchannel and permanent magnets
% %Plot of X-Y plane (top view) of the microchannel and permanent magnets
% figure_system_XY = figure('Name','System position X-Y plane','Color','white','Units', 'pixels','Position', [100 100 500 375]);
% hold on
% rectangle('Position',[-l_c -w_c 2*l_c 2*w_c], 'FaceColor',[1 0 0]); % microchannel topview, face color red
% rectangle('Position',[Magnet_1_Center(1)-Magnet_1_Dimension(1)/2 Magnet_1_Center(2)-Magnet_1_Dimension(2)/2 Magnet_1_Dimension(1) Magnet_1_Dimension(2)], 'EdgeColor',[0 0 0]); % Magnet 1 topview, face color black
%
% xlabel('X (m)')
% ylabel('Y (m)')
% set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',20,'FontWeight','Bold', 'LineWidth', 2, 'Box', 'on', 'XMinorTick','on','YMinorTick','on');
% axis equal
% set(gcf, 'PaperPositionMode', 'auto');
%
% %
% %Plot of Y-Z plane (cross-section view of the microchannel) and permanent magnets
% figure_system_YZ = figure('Name','System position Y-Z plane','Color','white','Units', 'pixels','Position', [100 100 500 375]);
% hold on
% rectangle('Position',[-w_c -t_c 2*w_c 2*t_c], 'FaceColor',[1 0 0]); % microchannel cross section view, face color red
% rectangle('Position',[Magnet_1_Center(2)-Magnet_1_Dimension(2)/2 Magnet_1_Center(3)-Magnet_1_Dimension(3)/2 Magnet_1_Dimension(2) Magnet_1_Dimension(3)], 'EdgeColor',[0 0 0]); % Magnet 1 topview, face color black
%
%
% xlabel('Y (m)')
% ylabel('Z (m)')
% set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',20,'FontWeight','Bold', 'LineWidth', 2, 'Box', 'on', 'XMinorTick','on','YMinorTick','on');
% axis equal
% set(gcf, 'PaperPositionMode', 'auto');

%% Plot tapered microchannel and magnet geometry
% figure('Name','Tapered Channel Geometry','Color','white','Units','pixels','Position',[100 100 600 400]);
% hold on;
%
% % Geometry parameters
% x_in  = -w_c;            % inlet position (left)
% x_out =  w_c;            % outlet position (right)
% w_in  = l_c/2;         % inlet half-width
% w_out =  l_out/2;         % outlet half-width
%
% % Define polygon points for tapered channel (top view)
% Xpoly = [-w_in -w_out  w_out  w_in];   % former width → X axis
% Ypoly = [ x_in  x_out   x_out  x_in];  % former length → Y axis
%
% fill(Xpoly, Ypoly, [1 0.8 0.8], 'EdgeColor','r', 'LineWidth',1.5);
% xlabel('X (m)');
% ylabel('Y (m)');
% title('Top View of Tapered Microchannel');
% axis equal;
% set(gca,'FontName','Arial','FontSize',16,'FontWeight','Bold','Box','on','LineWidth',1.5);
%
% % Draw the permanent magnet outline
% rectangle('Position', ...
%     [Magnet_1_Center(1)-Magnet_1_Dimension(1)/2, ...
%      Magnet_1_Center(2)-Magnet_1_Dimension(2)/2, ...
%      Magnet_1_Dimension(1), Magnet_1_Dimension(2)], ...
%      'EdgeColor',[0 0 0], 'LineWidth',1.5);
%
% legend('Tapered microchannel','Permanent magnet');
% hold off;


%% Solver for cell trajectories in microchannel
% Preallocating matrices for cell trajectories in X Y Z directions
XI=0;
Z_nb=1;Y_nb=1; % number of points on Y, Z direction
ZI = linspace(0,0,1); % Starting Z, unit m
YI=linspace(-w_c,-w_c,1);

for n_beads = 1:length(Par_beads)
    Par_bead = Par_beads{n_beads};
    Dia_cell=Par_bead(1);
    Vol_cell=Par_bead(2);
    Phi_p_cell = Par_bead(3);
    M_sc = Phi_p_cell*M_db;
    beta=1/3/pi/eta_f/Dia_cell; % Cell drag coefficienct

    for start_point_Y = 1:length(YI)
        YI_current = YI(start_point_Y);
        for start_point_Z = 1:length(ZI)
            ZI_current = ZI(start_point_Z);

            X0 = [XI YI_current ZI_current]; % Initial position (XI YI ZI) of cells in the stage II, unit m
            for cell_type = 1:length(Dia_cell)
                [t_cell,cell_traj] = ode45(@equations_motion,tspan,X0,options); % Solving motion of equation for cells
                t_cell_type1{start_point_Y,start_point_Z,cell_type} = t_cell;
                traj_cell_type{start_point_Y,start_point_Z,cell_type} = cell_traj;

                X_end1=cell_traj(end,1); %cell position near outlet
                Y_end1=cell_traj(end,2)*1e6 %cell position near outlet
                Z_end1=cell_traj(end,3); %cell position near outlet

            end
        end
    end
end

    %% Plot for cell trajectories in microchannel
    % % plot cell trajectories in X-Y plane of microchannel

     figure(19)
    for start_point_Y = 1:1:length(YI)
        for start_point_Z = 1:1:length(ZI)
            for cell_type = 1:length(Dia_cell)
                 Str(n_beads,cell_type) = cellstr(['D = ' num2str(Dia_cell)]);
                 H1(n_beads,cell_type) = plot(traj_cell_type{start_point_Y,start_point_Z,cell_type}(:,1),traj_cell_type{start_point_Y,start_point_Z,cell_type}(:,2),'LineWidth',2,'MarkerSize',10,'DisplayName',['beads #(' num2str(n_beads) ')']);

                 legend(H1,Str,'Location','northwest');
                hold on;
            end
        end
    end


    grid off;
    xlabel('X(m)');
    ylabel('Y(m)');
    xlim([-l_c l_c]);
    ylim([-w_c w_c]);
    %
    % set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 2, 'Box', 'on', 'XMinorTick','off','YMinorTick','off');
    % set(gcf, 'PaperPositionMode', 'auto');



%%plot cell trajectories in X-Z plane of microchannel
% figure(20)
% for start_point_Y = 1:1:length(YI)
%     for start_point_Z = 1:1:length(ZI)
%         for cell_type = 1:length(Dia_cell)
%                         Str(cell_type) = cellstr(['D = ' num2str(Dia_cell)]);
%             H1(cell_type) = plot(traj_cell_type{start_point_Y,start_point_Z,cell_type}(:,1),traj_cell_type{start_point_Y,start_point_Z,cell_type}(:,3),'r-','LineWidth',2,'MarkerSize',10);
%
%             hold on;
%             legend(H1,Str,'location','northwest');
%         end
%     end
% end
% grid off;
% xlabel('X(m)');
% ylabel('Z(m)');
% xlim([-l_c l_c]);
% ylim([-t_c t_c]);
% set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 2, 'Box', 'on', 'XMinorTick','off','YMinorTick','off');
% set(gcf, 'PaperPositionMode', 'auto');
%
% % % %
% % plot cell trajectories in Y-Z plane of microchannel
% figure(21)
% for start_point_Y = 1:1:length(YI)
%     for start_point_Z = 1:1:length(ZI)
%         for cell_type = 1:length(Dia_cell)
%                         Str(cell_type) = cellstr(['D = ' num2str(Dia_cell)]);
%             H1(cell_type) = plot(traj_cell_type{start_point_Y,start_point_Z,cell_type}(:,2),traj_cell_type{start_point_Y,start_point_Z,cell_type}(:,3),'r-','LineWidth',2,'MarkerSize',10);
%             hold on;
%             legend(H1,Str,'location','northwest');
%         end
%     end
% end
% grid off;
% xlabel('Y(m)');
% ylabel('Z(m)');
% xlim([-w_c w_c]);
% ylim([-t_c t_c]);
% set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',10,'FontWeight','Bold', 'LineWidth', 2, 'Box', 'on', 'XMinorTick','off','YMinorTick','off');
% set(gcf, 'PaperPositionMode', 'auto');
% print -depsc2 BeadTrajectories_Y_Z_plane.eps
%

% figure(22)
% for start_point_Y = 1:1:length(YI)
%     for start_point_Z = 1:1:length(ZI)
%         for cell_type = 1:length(Dia_cell)
%             H1(cell_type)=plot3(traj_cell_type{start_point_Y,start_point_Z,cell_type}(:,1),traj_cell_type{start_point_Y,start_point_Z,cell_type}(:,2),traj_cell_type{start_point_Y,start_point_Z,cell_type}(:,3),'r-','LineWidth',2,'MarkerSize',10);
%             hold on;
%             Str(cell_type) = cellstr(['D = ' num2str(Dia_cell)]);
%             legend(H1,Str,'location','northwest');
%         end
%     end
% end
% xlabel('X(m)');
% ylabel('Y(m)');
% zlabel('Z(m)');
% xlim([-l_c l_c]);
% ylim([-w_c w_c]);
% zlim([-t_c t_c]);
% set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',20,'FontWeight','Bold', 'LineWidth', 2, 'Box', 'on', 'XMinorTick','off','YMinorTick','off');
% set(gcf, 'PaperPositionMode', 'auto');



%%
% Function of motion equation
function motion = equations_motion(t,X)
global beta Dia_cell t_c w_c v_ave F_head_y l_c w_c t_c w_out l_out w_in

%% ------------------------------------------------------------------------
% 1. Geometry definition
% -------------------------------------------------------------------------
% The cell travels along the Y-direction.
% The channel width (along X) linearly tapers from the inlet to the outlet.

y_in  = -w_c;        % Bottom inlet position (Y = -w_c)
y_out =  w_c;        % Top outlet position  (Y = +w_c)
Y = X(2);            % Current Y position of the cell (in meters)

% Normalized position factor s in [0,1]
s = (Y - y_in) / (y_out - y_in);
s = max(0, min(1, s));   % Clamp to the valid range

% Compute local half-width of the channel in the X-direction
% At inlet: l_c_local = l_c (wide)
% At outlet: l_c_local = l_out (narrow)
l_c_local = (1 - s) * l_c + s * l_out;



%% Magnetic Force
[F_x, F_y, F_z, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = magnetic_field_force(X);

% Drag Coefficient (a particle near one wall, Brenner book, Lin et al PRE paper)
dis=abs((t_c-Dia_cell/2-abs(X(3))));
f_D_x_y=1/(1-(9/16)*(Dia_cell/(Dia_cell+2*dis))+(1/8)*(Dia_cell/(Dia_cell+2*dis))^3-(45/256)*(Dia_cell/(Dia_cell+2*dis))^4-(1/16)*(Dia_cell/(Dia_cell+2*dis))^5);
f_D_z=1/(1-(9/8)*(Dia_cell/(Dia_cell+2*dis))+(1/2)*(Dia_cell/(Dia_cell+2*dis))^3);
% f_D_x_y=1;
% f_D_z=1;
gamma_x_y=beta/f_D_x_y;
gamma_z=beta/f_D_z;


%% Velocity profile in X direction, based on Brody's equation
for i = 1:10
    n = i-1;
    R1 = (2*n+1)*pi/(2*t_c);
    R2 = (2*n+1)*pi/(4*t_c);
    b(i) = 1/(2*n+1)^4*(1-tanh(R2*2*l_c_local)/(R2*2*l_c_local));
    a(i) = (-1)^n/(2*n+1)^3*(1-cosh(R1*X(2))/cosh(R2*2*l_c_local))*cos(R1*(X(3)));

end
v_f_y_z = pi*v_ave/2 * sum(a)/sum(b);

%% Cell stiffness force (Hertz–Tatara–Hyperelastic)
%% Cell stiffness force (continuous soft→hard contact)
[F_stiff, modelID] = contact_force_auto_tapered(l_c_local, t_c);

% stiffness force in the opsite direction of the external forces
driveY = F_head_y + F_y;
sgnY   = sign(driveY);  if sgnY == 0, sgnY = 1; end

F_stiff_y = -sgnY * F_stiff;   % move only in Y direction
F_stiff_z = 0;                 % neglect z direction



%% Total forces
F_total_y = F_y + F_head_y + F_stiff_y;
F_total_z = F_z - abs(F_stiff_z);


%% Set up boundary conditions
% --- Y direction (lateral) ---

motion_y = F_total_y * gamma_x_y;


% prevent reverse motion
if F_total_y < 0
    motion_y = 0;
end


% --- Z direction (vertical) ---
at_wall_z = abs(X(3)) >= (t_c - Dia_cell/2);

if at_wall_z && (sign(F_total_z) * sign(X(3))) > 0
    motion_z = 0;
else
    motion_z = F_total_z * gamma_z;
end

if F_total_z < 0
    motion_z = 0;
end


% Equations of motion
motion = [0; motion_y; 0];
end



%% Function of magnetic field and force
function [F_x, F_y, F_z, H_x, H_y, H_z, H_xx, H_xy, H_xz, H_yx, H_yy, H_yz, H_zx, H_zy, H_zz] = magnetic_field_force(X)
global Par_magnets u_0 M_sf Vol_cell d_pf ln_sigma_pf M_d ka T d_pb ln_sigma_pb M_db M_sc

dX = 1e-6;

Hx_before_summation = zeros(length(Par_magnets),7);
Hy_before_summation = zeros(length(Par_magnets),7);
Hz_before_summation = zeros(length(Par_magnets),7);


for n_magnets = 1:length(Par_magnets)
    Par_magnet = Par_magnets{n_magnets};
    M_s = Par_magnet(1,2)/u_0;
    l = Par_magnet(2,1)/2;
    w = Par_magnet(2,2)/2;
    h = Par_magnet(2,3)/2;
    x0 = Par_magnet(3,1);
    y0 = Par_magnet(3,2);
    z0 = Par_magnet(3,3);
    rot_x = Par_magnet(4,1);
    rot_y = Par_magnet(4,2);
    rot_z = Par_magnet(4,3);

    Xm=[-l, l];
    Ym=[-w, w];
    Zm=[-h, h];

    Mx = [1 0 0; 0 cosd(rot_x) -sind(rot_x); 0 sind(rot_x) cosd(rot_x)];
    My = [cosd(rot_y) 0 sind(rot_y); 0 1 0; -sind(rot_y) 0 cosd(rot_y)];
    Mz = [cosd(rot_z) -sind(rot_z) 0; sind(rot_z) cosd(rot_z) 0; 0 0 1];

    Mxyz = Mz*My*Mx;
    Mzyx = inv(Mxyz);

    for m_gradient_loop = 1:7
        XYZ = X;
        mm = mod(m_gradient_loop,7);
        if mm>0
            p = ceil(mm/2);
            q = mod(mm,2);
            XYZ(p) = X(p)+(-1)^q*dX;
        else
            XYZ = X;
        end

        xyz = Mxyz*(XYZ'-[x0 y0 z0])';
        %xyz = Mxyz*(XYZ'-[x0 y0 z0])';
        x = xyz(1);
        y = xyz(2);
        z = xyz(3);

        for k=1:2
            for m=1:2
                h_x(k,m)=(-1)^(k+m)*log(((z-(-h))+((x-Xm(m))^2+(z-(-h))^2+(y-Ym(k))^2)^0.5)/((z-(+h))+((x-Xm(m))^2+(z-(+h))^2+(y-Ym(k))^2)^0.5));
            end
        end

        for k=1:2
            for m=1:2
                h_z(k,m)=(-1)^(k+m)*log(((x-(-l))+((y-Ym(k))^2+(x-(-l))^2+(z-Zm(m))^2)^0.5)/((x-(+l))+((y-Ym(k))^2+(x-(+l))^2+(z-Zm(m))^2)^0.5));
            end
        end

        for k=1:2
            for n=1:2
                for m=1:2
                    h_y(k,n,m)=(-1)^(k+n+m)*atan((x-Xm(n))*(z-Zm(m))/((y-Ym(k))*((x-Xm(n))^2+(y-Ym(k))^2+(z-Zm(m))^2)^0.5));
                end
            end
        end
        hx = M_s/(4*pi)*nansum(nansum(h_x));
        hz = M_s/(4*pi)*nansum(nansum(h_z));
        hy = M_s/(4*pi)*nansum(nansum(nansum(h_y)));
        H_tran = Mzyx*[hx hy hz]';
        Hx_before_summation(n_magnets,m_gradient_loop) = H_tran(1);
        Hy_before_summation(n_magnets,m_gradient_loop) = H_tran(2);
        Hz_before_summation(n_magnets,m_gradient_loop) = H_tran(3);
    end
end

% field gradient calculation (central differencing method)
H_xx = (sum(Hx_before_summation(:,2))-sum(Hx_before_summation(:,1)))/2/dX;
H_xy = (sum(Hx_before_summation(:,4))-sum(Hx_before_summation(:,3)))/2/dX;
H_xz = (sum(Hx_before_summation(:,6))-sum(Hx_before_summation(:,5)))/2/dX;
H_yx = (sum(Hy_before_summation(:,2))-sum(Hy_before_summation(:,1)))/2/dX;
H_yy = (sum(Hy_before_summation(:,4))-sum(Hy_before_summation(:,3)))/2/dX;
H_yz = (sum(Hy_before_summation(:,6))-sum(Hy_before_summation(:,5)))/2/dX;
H_zx = (sum(Hz_before_summation(:,2))-sum(Hz_before_summation(:,1)))/2/dX;
H_zy = (sum(Hz_before_summation(:,4))-sum(Hz_before_summation(:,3)))/2/dX;
H_zz = (sum(Hz_before_summation(:,6))-sum(Hz_before_summation(:,5)))/2/dX;
% field summation of multiple magnets
H_x = sum(Hx_before_summation(:,7));
H_y = sum(Hy_before_summation(:,7));
H_z = sum(Hz_before_summation(:,7));
H_ex=norm([H_x H_y H_z]);


%% Magnetic force without considering demagnetization
% ferrofluid magnetization considering just a mean particle diameter
%Mf = M_sf*(coth(alpha_f*H_ex)-1/alpha_f/H_ex);

%ferrofluid magnetization considering a log-normal distribution of nanoparticles, when H_ex is zero, magnetization is set to zero, because coth(0)-1/0 is not defined in MATLAB.
if (H_ex == 0)
    Mf = 0;
else
    Mf = M_sf*(integral(@(d_m) ((1./(sqrt(2.*pi).*d_m.*ln_sigma_pf)).*exp(-((log(d_m./d_pf)).^2./(2.*ln_sigma_pf.^2)))).*(coth(pi.*u_0.*M_d.*d_m.^3.*H_ex./(6.*ka.*T)) - (6.*ka.*T)./(pi.*u_0.*M_d.*d_m.^3.*H_ex)), 0, Inf));
end
% microbead magnetization considering just a mean particle diameter
%Mp = M_sc*(coth(alpha_p*H_ex)-1/alpha_p/H_ex);

%magnetic bead magnetization considering a log-normal distribution of nanoparticles, when H_ex is zero, magnetization is set to zero, because coth(0)-1/0 is not defined in MATLAB.
if (H_ex == 0)
    Mp =0;
else
    Mp = M_sc*(integral(@(d_m) ((1./(sqrt(2.*pi).*d_m.*ln_sigma_pb)).*exp(-((log(d_m./d_pb)).^2./(2.*ln_sigma_pb.^2)))).*(coth(pi.*u_0.*M_db.*d_m.^3.*H_ex./(6.*ka.*T)) - (6.*ka.*T)./(pi.*u_0.*M_db.*d_m.^3.*H_ex)), 0, Inf));
end

% difference between ferrofluid and bead magnetization
M = -Mf;

if (H_ex == 0)
    M_x=0;
    M_y=0;
    M_z=0;
else
    M_x=M*H_x/H_ex;
    M_y=M*H_y/H_ex;
    M_z=M*H_z/H_ex;
end

F_x=u_0*Vol_cell*(M_x*H_xx+M_y*H_xy+M_z*H_xz);
F_y=u_0*Vol_cell*(M_x*H_yx+M_y*H_yy+M_z*H_yz);
F_z=u_0*Vol_cell*(M_x*H_zx+M_y*H_zy+M_z*H_zz);


end

%% Function of end of channel detection

function [value,isterminal,direction] = events_end_of_channel(t,X,p)
global XE
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.

value = X(2)-XE;
isterminal = [1];   % stop the integration
direction = [0];   % negative direction
end


%% Stiffness force with tapered geometry
function [Fstiff, modelID] = contact_force_auto_tapered(l_c_local, t_c)
% Computes total stiffness/contact force (Y and Z walls)
% for a channel whose width tapers along X.

global Dia_cell E_cell nu

R = Dia_cell / 2;                  % Cell radius
Estar = E_cell / (1 - nu^2);       % Effective modulus for nearly incompressible material

    function [F, id] = one_dir_force(Wdir)
        % ---------------------------------------------------------------
        % Compute wall-contact stiffness in one direction
        % Input:  Wdir = local half-gap between walls (m)
        % Output: F (N), id (model regime string)
        % ---------------------------------------------------------------

        % --- Physical overlap distance ---
        % Positive delta means compression, negative means free gap
        margin = 0.2e-6;                       % 0.2 µm soft-contact margin
        delta = R - (Wdir + margin);           % Overlap distance

        % --- Free region: no contact ---
        if delta <= 0
            F = 0;
            id = "no-contact";
            return;
        end

        % --- Hertzian contact regime ---
        x = delta / Dia_cell;                  % Normalized deformation
        a = 0.5 * sqrt( max(Dia_cell^2 - (2*Wdir)^2, 0) );
        Fh = (4/3) * Estar * sqrt(R) * delta^(3/2);
        if x < 0.10
            F = Fh;
            id = "Hertz";
            return;
        end

        % --- Tatara regime (moderate deformation) ---
        fa = (1 - 1.5*(a/Dia_cell) + 0.5*(a/Dia_cell)^3);
        Ft = Fh / max(fa, eps);
        if x <= 0.30
            F = Ft;
            id = "Tatara";
            return;
        end

        % --- Hyperelastic–Tatara regime (large deformation) ---
        A = (1 - x)^2 / (1 - x + x^2)^3;
        B = (1 - x/3) / (1 - x + x^2)^3;
        term1 = ((3*(1 - nu^2)*A)/2) * (delta/max(a,eps)) * (1 + 2*B*(a^2/Dia_cell^2));
        term2 = 2*A*sqrt(pi)*delta*(1 + 4*B*a^2/(5*Dia_cell^2)) * ...
            (1 - 1.5*(a/Dia_cell) + 0.5*(a/Dia_cell)^3);
        C_hyper = term1 - term2;
        F = E_cell / max(C_hyper, eps);
        id = "Hyperelastic–Tatara";
    end

% Compute forces for lateral (tapered) and vertical walls
[Fy, id_y] = one_dir_force(l_c_local);   % Lateral Y-direction (channel taper)
[Fz, id_z] = one_dir_force(t_c);         % Top–bottom Z-direction

% Total contact force (sum of components)
Fstiff = Fy + Fz;
modelID.y = id_y;
modelID.z = id_z;

end