%% FILE PATH
path_to_folder_with_wax_data_files = 'C:\Users\la3314de\Downloads\WAXS_data\';
[I_bsd, I_bck, Im, Im0] = load_wax_data(path_to_folder_with_wax_data_files);

%I_bsd - The beam-stop diode value during scattering.
%I_bck - The beam-stop diode value during background collection.
%Im - The detector images during scattering.
%Im0 - The detector images during background collection.

%% CHECK NOT NORMALIZED PICS
 figure(1)
for i=1:6
axis = subplot(2,3,i);
 pcolor( axis, Im{i} )
 shading interp
 daspect([1 1 1])
 title(['Load step ',num2str(i)])
end

%% FUNCTION TO READ IN load_wax_data.m

function [diod_values, background_diod_values, scattering_images, background_images] = load_wax_data(path_to_folder_with_wax_data_files)
% ------------------------------------------------------------------- -------------
% The output is:
% diod_values - The beamstop diod value during scattering.
% background_diod_values - The beamstop diod value during background  collection.
% scattering_images - The detector images during scattering.
% background_images - The detector images during background collection.
% --------------------------------------------------------------------------------
diod_values = {};
background_diod_values = {};
scattering_images = {};
background_images = {};
stem = [path_to_folder_with_wax_data_files, 'pc_05s_03_130604_'];
names = [1,9,12,14,16,23];
for i=1:6
 
 I = double(imread(join([stem, num2str(names(i),'%03.f'),'_007.tif'],"")) );
 fileID = fopen(join([stem,num2str(names(i),'%03.f'),'_007.tif.hdr'],""),'r');
 s1 = fscanf(fileID, '%s');
 s2 = strsplit(s1,'BSd=');
 s3 = strsplit( s2{2}, 'Vacuum');
 diodI = str2double( s3{1} );
 fclose(fileID);
 
 I0 = double(imread(join([stem,num2str(names(i),'%03.f'),'_000.tif'],"")) );
 fileID = fopen(join([stem,num2str(names(i),'%03.f'),'_000.tif.hdr'],""),'r');
 s1 = fscanf(fileID, '%s');
 s2 = strsplit(s1,'BSd=');
 s3 = strsplit( s2{2}, 'Vacuum');
 diodI0 = str2double( s3{1} );
 fclose(fileID);

 diod_values{end+1} = diodI;
 background_diod_values{end+1} = diodI0;
 scattering_images{end+1} = I;
 background_images{end+1} = I0;
end

end

%% READ IN THE POINTS
xpoints = readmatrix('WAXS_data/X-ray_points.txt');
force_disp = readmatrix('WAXS_data/load_data.txt');

%% CLEAN AND NORMALIZE THE IMAGES

Im_corr = {}; % Create an empty list to store corrected images
global_min = inf;  % Start with a very high number
global_max = -inf; % Start with a very low number

for i = 1:6
    % Get the Raw Scattering Image and its Diode Value
    Scatter_Image = Im{i};       
    Scatter_Diode = I_bsd{i};    
    
    % Get the Background Image and its Diode Value
    Backgr_Image = Im0{i};       
    Backgr_Diode = I_bck{i};     
    
    % Apply the Formula: (Im / I_bsd) - (Im0 / I_bck)
    % Note: We divide the whole image matrix by the single diode scalar
    Corrected_Image = (Scatter_Image / Scatter_Diode) - (Backgr_Image / Backgr_Diode);
    
    % 4. Store the result
    Im_corr{i} = Corrected_Image;

    % -- UPDATE GLOBAL LIMITS --
    % Check if this image has a new max or min
    current_max = max(Corrected_Image(:));
    current_min = min(Corrected_Image(:));
    
    if current_max > global_max
        global_max = current_max;
    end
    if current_min < global_min
        global_min = current_min;
    end
end

figure(2); clf; 
for i = 1:6
    ax = subplot(2, 3, i);
    pcolor(ax, Im_corr{i})      
    shading interp
    daspect([1 1 1])
    title(['Step ', num2str(i)])
    
    % -- APPLY THE GLOBAL LIMITS --
    % This forces every subplot to use the same scale
    clim([global_min, global_max]); 
    
end

% Add a single colorbar to the side 
hp4 = get(subplot(2,3,6),'Position'); 
colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.02  hp4(4)*2.35])

% Save the 1D Scattering Patterns (Figure 2)
exportgraphics(figure(2), 'WAXS_1D_Patterns.png', 'Resolution', 300);

%% --- PARAMETER SETUP  ---
lambda = 0.91e-10;      % Wavelength in meters (0.91 Angstrom)
D = 0.3768;               % Sample-to-detector distance in meters (WAXS)
pixel_size = 172e-6;    % Pixel size in meters (Standard Pilatus detector)

% Get image size
[rows, cols] = size(Im_corr{1});

% Define Beam Center
cen_y = 575.54; %This is in px. 
cen_x = 329.45; 

%% --- GENERATE Q-MAP ---
% We calculate the 'q' value for every pixel in the image.
% q = (4*pi/lambda) * sin(theta)
% tan(2*theta) = radius / distance

[X, Y] = meshgrid(1:cols, 1:rows);
R_px = sqrt((X - cen_x).^2 + (Y - cen_y).^2); % Distance from center in pixels
R_m = R_px * pixel_size;                      % Distance from center in meters

theta = 0.5 * atan(R_m / D);                  % Half scattering angle (radians)
Q_map = (4 * pi / lambda) * sin(theta);       % Scattering vector map (1/m)
Q_nm = Q_map / 1e9;                           % Convert to 1/nm for easier reading

%% --- EXTRACT 1D PROFILES USING INTERP2 ---
% We create the coordinate grid once, as it is the same for all images.
[XI, YI] = meshgrid(1:cols, 1:rows);

% Define the lines we want to extract (Radial lines from center to edge)
% We calculate 'r' (radius in pixels) to generate coordinates.

% -- Horizontal Line (Perpendicular to Load) --
% Going from center to the right edge
len_horz = cols - ceil(cen_x);       % Max length to the right
r_vec_horz = 0:len_horz;             % Radius vector (pixels)
x_line_horz = cen_x + r_vec_horz;    % X coordinates
y_line_horz = ones(size(x_line_horz)) * cen_y; % Y coordinates (constant)

% -- Vertical Line (Parallel to Load) --
% Going from center to the top (or bottom, depending on where the data is good)
% Let's go 'up' (decreasing Y) as an example. If blocked, change sign to +
len_vert = ceil(cen_y);              % Max length to top
r_vec_vert = 0:len_vert;             % Radius vector (pixels)
x_line_vert = ones(size(r_vec_vert)) * cen_x;  % X coordinates (constant)
y_line_vert = cen_y - r_vec_vert;    % Y coordinates (going up)

% --- Convert Pixel Radius to q for the x-axis ---
% We can convert the r_vectors directly to q.
% Function: q = (4pi/lambda) * sin( 0.5 * atan( r*pixel_size / D ) )

% Horizontal q-axis
R_m_horz = r_vec_horz * pixel_size;
theta_horz = 0.5 * atan(R_m_horz / D);
q_axis_horz = (4 * pi / lambda) * sin(theta_horz) / 1e9; % /1e9 for nm^-1

% Vertical q-axis
R_m_vert = r_vec_vert * pixel_size;
theta_vert = 0.5 * atan(R_m_vert / D);
q_axis_vert = (4 * pi / lambda) * sin(theta_vert) / 1e9;


%% --- 1D INTENSITIES VS Q ---
figure(3); clf;
% Make the figure wider: [Left, Bottom, Width, Height]
set(gcf, 'Position', [100, 100, 1600, 550]); 

% Setup colors for the 6 load steps
colors = jet(6); 

% -- PREPARE SUBPLOT 1 (Vertical) --
subplot(1,2,1); hold on; box on;
title('Vertical Cut (Parallel to Load)');
xlabel('q (nm^{-1})'); ylabel('Intensity (a.u.)');
xlim([-2, 25]);   % Focus on relevant q-range
ylim([0, 0.35]); % Force Y-axis to start at 0 (Hides negative noise)

% -- PREPARE SUBPLOT 2 (Horizontal) --
subplot(1,2,2); hold on; box on;
title('Horizontal Cut (Perpendicular)');
xlabel('q (nm^{-1})'); ylabel('Intensity (a.u.)');
xlim([-2, 25]);   
ylim([0, 0.35]); % Force Y-axis to start at 0

for i = 1:6
    Current_Image = Im_corr{i};
    
    % Interpolate Horizontal Intensity
    I_profile_horz = interp2(XI, YI, Current_Image, x_line_horz, y_line_horz);
    
    % Interpolate Vertical Intensity
    I_profile_vert = interp2(XI, YI, Current_Image, x_line_vert, y_line_vert);
    
    % Plot Vertical (Subplot 1)
    subplot(1,2,1);
    plot(q_axis_vert, I_profile_vert, 'Color', colors(i,:), 'LineWidth', 1.5, ...
        'DisplayName', ['Step ' num2str(i)]);
    
    % Plot Horizontal (Subplot 2)
    subplot(1,2,2);
    plot(q_axis_horz, I_profile_horz, 'Color', colors(i,:), 'LineWidth', 1.5, ...
        'DisplayName', ['Step ' num2str(i)]);
end

% Add Legends
subplot(1,2,1); legend show;
subplot(1,2,2); legend show;

% Save the 1D Scattering Patterns (Figure 2)
exportgraphics(figure(3), 'WAXS_q_vs_intensity.png', 'Resolution', 300);

%% --- AZIMUTHAL POLAR PLOTS (Circular Intensity Profile) ---

% The distance from the center represents Intensity.
% Question 7: Look at variation around a circle at constant q.
% We need to find the q-radius of the main peak first.

% Define Parameters
target_q = 6.085;  % PEAK POSITION 
q_tol = 0.02;       % Width of the integration ring (0.04 for Peak, 0.4 Peak 2)
colors = jet(6);   % Color scheme

figure(4); clf;
% Create a polar axes
pax = polaraxes; 
hold on;
title(['Azimuthal Intensity Distribution at q \approx ' num2str(target_q) ' nm^{-1}']);

peak_anisotropy = zeros(1,6); 

for i = 1:6
    curr_Im = Im_corr{i};
    
    % --- Extract the Ring ---
    ring_mask = (Q_nm >= (target_q - q_tol)) & (Q_nm <= (target_q + q_tol));
    
    % Get coordinates and intensity of pixels in the ring
    [y_ring, x_ring] = find(ring_mask);
    I_ring_vals = curr_Im(ring_mask);
    
    % --- Calculate Angles (Radians) ---
    % atan2 returns values between -pi and +pi
    angles_rad = atan2(y_ring - cen_y, x_ring - cen_x);
    
    % --- Sort and Smooth ---
    [sorted_angles, sort_idx] = sort(angles_rad);
    sorted_I = I_ring_vals(sort_idx);
    
    % Smooth data to get a clean "shape" instead of noisy spikes
    smooth_I = smoothdata(sorted_I, 'gaussian', 100); 
    
    % --- Close the Loop (Connect last point to first) ---
    % We append the first point to the end of the arrays
    plot_angles = [sorted_angles; sorted_angles(1)];
    plot_intensity = [smooth_I; smooth_I(1)];
    
    % --- Polar Plot ---
    % polarplot(theta_radians, rho_intensity)
    polarplot(pax, plot_angles, plot_intensity, 'Color', colors(i,:), ...
        'LineWidth', 2, 'DisplayName', ['Step ' num2str(i)]);
    
    % --- Calculate Anisotropy Factor ---
    %peak_anisotropy(i) = max(smooth_I) - min(smooth_I);
end

legend('Location', 'northeastoutside'); 

% Optional: Adjust limits to make the shape clearer
rlim([0, inf]);

% Save Figure
exportgraphics(figure(4), 'WAXS_azimutha6.png', 'Resolution', 300);

%% ---  MACROSCOPIC CURVE ---
%  Force-Displacement curve.

figure(5); clf;

% -- LEFT AXIS: MACROSCOPIC DATA --
yyaxis left
% Plot the full continuous Load Curve (Black Line)
% Column 1 = Displacement, Column 2 = Force
plot(force_disp(:,1), force_disp(:,2), 'k-', 'LineWidth', 1);
ylabel('Force (N)');
xlabel('Displacement (mm)');
hold on;

% Plot the 6 WAXS measurement points (Blue Dots)
% We plot Column 1 (Disp) vs Column 2 (Force) directly.
waxs_disp = xpoints(:,1); 
waxs_force = xpoints(:,2);

plot(waxs_disp, waxs_force, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);

% -- RIGHT AXIS: MICROSCOPIC ANISOTROPY --
%yyaxis right
% Plot the Anisotropy factor calculated in the loop
% X-axis: The displacement at the 6 points
% Y-axis: The calculated anisotropy
%plot(waxs_disp, peak_anisotropy, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');

%ylabel('WAXS Anisotropy (Int_{max} - Int_{min})');
ax = gca;
ax.YColor = 'k'; % Make right axis red to match the curve

%title('Macroscopic Load vs Microscopic Anisotropy');
legend('Force Curve', 'WAXS Steps', 'Location', 'best'); %'Anisotropy Evolution',
grid on;

exportgraphics(figure(5), 'WAXS_Macro_vs_Micro.png', 'Resolution', 300);

%% --- PEAK AMPLITUDE VS AZIMUTH (Cartesian Plot) ---
% Peak Intensity vs Angle.

figure(6); clf;
set(gcf, 'Position', [100, 100, 1000, 600]); % Nice wide window
hold on; box on;
%title('Anisotropy Evolution: Peak Amplitude vs. Azimuth');
xlabel('Azimuthal Angle (Degrees)');
ylabel('Peak Intensity (a.u.)');
xlim([-180 180]); % -180 is Left, 0 is Right, 90 is Top
xticks(-180:45:180); % Ticks every 45 degrees
grid on;

% --- Parameters ---
n_sectors = 72;       % Divide 360 degrees into 72 steps (5-degree resolution)
q_search_min = 10;    % Look for the peak inside this range (Amorphous Halo)
q_search_max = 20;    

% Generate angle bins (-pi to +pi)
edges = linspace(-pi, pi, n_sectors + 1);
sector_centers = (edges(1:end-1) + edges(2:end)) / 2;
sector_deg = rad2deg(sector_centers); % Convert to degrees for plotting

colors = jet(6);

% --- Pre-calculate Angle Map ---
% We calculate the angle for every pixel once to save time
[X_coords, Y_coords] = meshgrid(1:cols, 1:rows);
Angle_Map = atan2(Y_coords - cen_y, X_coords - cen_x); % Returns -pi to +pi

for i = 1:6
    curr_Im = Im_corr{i};
    amplitudes = zeros(1, n_sectors);
    
    for k = 1:n_sectors
        % 1. Create a Mask for this specific "Pizza Slice" (Wedge)
        % Logic: Pixel angle is between edge A and edge B, AND Pixel is in the peak q-range
        wedge_mask = (Angle_Map >= edges(k)) & (Angle_Map < edges(k+1)) & ...
                     (Q_nm > q_search_min) & (Q_nm < q_search_max);
        
        % 2. Extract pixels
        vals = curr_Im(wedge_mask);
        
        % 3. Find the Peak Amplitude in this wedge
        if isempty(vals)
            amplitudes(k) = NaN; % Handle empty corners
        else
            % We take the max value. 
            amplitudes(k) = max(vals); 
        end
    end
    
    % 4. Fix "Dips" (Detector Gaps) using Smoothing
    % If a wedge hit a detector gap, the amplitude drops. We smooth over angles to fix this.
    amplitudes_smooth = smoothdata(amplitudes, 'gaussian', 6);
    
    % 5. Plot
    plot(sector_deg, amplitudes_smooth, 'Color', colors(i,:), ...
        'LineWidth', 2, 'DisplayName', ['Step ' num2str(i)]);
end

legend('Location', 'northeastoutside');
exportgraphics(figure(6), 'WAXS_Peak_Amplitude_vs_Azimuth.png', 'Resolution', 300);

%% --- PEAK POSITION (Q-VALUE) VS AZIMUTH ---

figure(7); clf;
set(gcf, 'Position', [100, 100, 1000, 600]); 
hold on; box on;
%title('Structural Distortion: Peak Position (q) vs. Azimuth');
xlabel('Azimuthal Angle (Degrees)');
ylabel('Peak Position q (nm^{-1})');
xlim([-180 180]); 
xticks(-180:45:180);
grid on;

% --- Parameters ---
n_sectors = 72;       % 5-degree steps
q_search_min = 10;    % Range to look for the Halo
q_search_max = 20;

% Angle Bins
edges = linspace(-pi, pi, n_sectors + 1);
sector_centers = (edges(1:end-1) + edges(2:end)) / 2;
sector_deg = rad2deg(sector_centers);

colors = jet(6);

% --- Pre-calculate Angle Map (if not already done) ---
[X_coords, Y_coords] = meshgrid(1:cols, 1:rows);
Angle_Map = atan2(Y_coords - cen_y, X_coords - cen_x); 

for i = 1:6
    curr_Im = Im_corr{i};
    peak_q_vals = zeros(1, n_sectors);
    
    for k = 1:n_sectors
        % 1. Create Mask for this wedge
        wedge_mask = (Angle_Map >= edges(k)) & (Angle_Map < edges(k+1)) & ...
                     (Q_nm > q_search_min) & (Q_nm < q_search_max);
        
        % 2. Extract q and Intensity for pixels in this wedge
        I_vals = curr_Im(wedge_mask);
        q_vals = Q_nm(wedge_mask);
        
        if isempty(I_vals)
            peak_q_vals(k) = NaN;
        else
            % 3. Find the q-value where Intensity is Maximum
            [~, max_idx] = max(I_vals);
            peak_q_vals(k) = q_vals(max_idx);
        end
    end
    
    % 4. Smooth the curve
    peak_q_smooth = smoothdata(peak_q_vals, 'gaussian', 10);
    
    % 5. Plot
    plot(sector_deg, peak_q_smooth, 'Color', colors(i,:), ...
        'LineWidth', 2, 'DisplayName', ['Step ' num2str(i)]);
end

legend('Location', 'northeastoutside');
exportgraphics(figure(7), 'WAXS_Peak_Position_vs_Azimuth.png', 'Resolution', 300);
