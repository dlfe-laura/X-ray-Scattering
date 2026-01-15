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
% Load WAX data for Modern exp. mech course LAB4. The output is:
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

figure(2); clf; % Clear figure to ensure clean slate
for i = 1:6
    ax = subplot(2, 3, i);
    pcolor(ax, Im_corr{i})      
    shading interp
    daspect([1 1 1])
    title(['Step ', num2str(i)])
    
    % -- APPLY THE GLOBAL LIMITS --
    % This forces every subplot to use the same scale
    clim([global_min, global_max]); 
    
    % Optional: If a few "hot pixels" make the image look black, 
    % try manually capping the max limit (uncomment below to try):
    % clim([0, global_max * 0.1]); 
end

% Add a single colorbar to the side (optional, looks cleaner)
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
% NOTE: In real experiments, this is calibrated
cen_y = 575.54; %This is in px.  rows / 2
cen_x = 329.45; %cols / 2

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

%% --- 3. EXTRACT 1D PROFILES USING INTERP2 ---
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



%% --- 4. PLOT THE INTENSITIES VS Q ---
figure(3); clf;
% Make the figure wider: [Left, Bottom, Width, Height]
set(gcf, 'Position', [100, 100, 1600, 550]); 

%sgtitle('1D Intensity Profiles (Interpolated)');

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

%% --- 5. PEAK ANALYSIS (Brief check) ---
% This section helps you answer Question 5 ("Are there any peaks?")
% We look for the max intensity in the amorphous halo region (e.g. q=10-20)

% Just taking the last load step as an example
mask_peak = q_axis_horz > 10 & q_axis_horz < 20;
[max_val, idx] = max(I_profile_horz(mask_peak));
q_peak_val = q_axis_horz(find(mask_peak, 1) + idx - 1);

fprintf('Step 6 Horizontal Peak Position: q = %.2f nm^-1\n', q_peak_val);
d_spacing = 2*pi / q_peak_val;
fprintf('Corresponds to physical spacing d = %.2f nm\n', d_spacing);

%% --- 4. EXTRACT, PLOT PROFILES AND FIND PEAKS ---
figure(3); clf;
set(gcf, 'Position', [100, 100, 1600, 550]); 

% Setup colors and legends
colors = jet(6); 
legend_list = cell(1,6);

% Prepare Subplots
ax1 = subplot(1,2,1); hold on; box on;
title('Vertical Cut (Parallel to Load)');
xlabel('q (nm^{-1})'); ylabel('Intensity (a.u.)');
xlim([-2, 25]); ylim([0, 0.35]);

ax2 = subplot(1,2,2); hold on; box on;
title('Horizontal Cut (Perpendicular)');
xlabel('q (nm^{-1})'); ylabel('Intensity (a.u.)');
xlim([-2, 25]); ylim([0, 0.35]);

% Initialize arrays to store peak data for the report
% Columns: [q_position, intensity, d_spacing]
results_vert = zeros(6, 3);
results_horz = zeros(6, 3);

% Define search range for the amorphous halo (e.g., 10 to 20 nm^-1)
q_search_min = 10;
q_search_max = 20;

for i = 1:6
    Current_Image = Im_corr{i};
    legend_list{i} = ['Step ' num2str(i)];
    
    % --- 1. EXTRACT DATA ---
    % Use median filter to prevent "dips" to zero caused by detector gaps
    I_profile_horz_raw = interp2(XI, YI, Current_Image, x_line_horz, y_line_horz);
    I_profile_horz = movmedian(I_profile_horz_raw, 5); % Cleans noise
    
    I_profile_vert_raw = interp2(XI, YI, Current_Image, x_line_vert, y_line_vert);
    I_profile_vert = movmedian(I_profile_vert_raw, 5); % Cleans noise
    
    % --- 2. FIND PEAKS ---
    % Vertical Peak
    mask_v = q_axis_vert > q_search_min & q_axis_vert < q_search_max;
    q_roi_v = q_axis_vert(mask_v);
    I_roi_v = I_profile_vert(mask_v);
    
    [max_I_v, idx_v] = max(smoothdata(I_roi_v, 'gaussian', 10)); % Smooth before finding max
    peak_q_v = q_roi_v(idx_v);
    
    % Horizontal Peak
    mask_h = q_axis_horz > q_search_min & q_axis_horz < q_search_max;
    q_roi_h = q_axis_horz(mask_h);
    I_roi_h = I_profile_horz(mask_h);
    
    [max_I_h, idx_h] = max(smoothdata(I_roi_h, 'gaussian', 10)); 
    peak_q_h = q_roi_h(idx_h);
    
    % Store Results
    results_vert(i,:) = [peak_q_v, max_I_v, 2*pi/peak_q_v];
    results_horz(i,:) = [peak_q_h, max_I_h, 2*pi/peak_q_h];
    
    % --- 3. PLOT TRACES AND MARKERS ---
    % Vertical
    plot(ax1, q_axis_vert, I_profile_vert, 'Color', colors(i,:), 'LineWidth', 1.5);
    plot(ax1, peak_q_v, max_I_v, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
    
    % Horizontal
    plot(ax2, q_axis_horz, I_profile_horz, 'Color', colors(i,:), 'LineWidth', 1.5);
    plot(ax2, peak_q_h, max_I_h, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
end

% Add Legends
legend(ax1, legend_list, 'Location', 'northeast');
legend(ax2, legend_list, 'Location', 'northeast');

% Save Figure
exportgraphics(figure(3), 'WAXS_q_vs_intensity_with_peaks.png', 'Resolution', 300);

%% --- 5. PRINT PEAK ANALYSIS TABLE ---
fprintf('\n=================================================================\n');
fprintf('                WAXS PEAK ANALYSIS REPORT                        \n');
fprintf('=================================================================\n');
fprintf('| Step | Dir  | Peak q (nm^-1) | d-spacing (nm) | Intensity (a.u) |\n');
fprintf('|------|------|----------------|----------------|-----------------|\n');

for i = 1:6
    % Print Vertical Data
    fprintf('|  %d   | Vert |     %6.2f     |     %6.3f      |     %6.3f      |\n', ...
        i, results_vert(i,1), results_vert(i,3), results_vert(i,2));
    
    % Print Horizontal Data
    fprintf('|  %d   | Horz |     %6.2f     |     %6.3f      |     %6.3f      |\n', ...
        i, results_horz(i,1), results_horz(i,3), results_horz(i,2));
    fprintf('|------|------|----------------|----------------|-----------------|\n');
end
fprintf('\n');

%% --- 4. EXTRACT, PROCESS AND FIND MULTIPLE PEAKS ---
figure(3); clf;
set(gcf, 'Position', [100, 100, 1600, 550]); 

% Setup colors and legends
colors = jet(6); 
legend_list = cell(1,6);

% Prepare Subplots
ax1 = subplot(1,2,1); hold on; box on;
title('Vertical Cut (Parallel to Load)');
xlabel('q (nm^{-1})'); ylabel('Intensity (a.u.)');
xlim([-2, 25]); ylim([0, 0.35]);

ax2 = subplot(1,2,2); hold on; box on;
title('Horizontal Cut (Perpendicular)');
xlabel('q (nm^{-1})'); ylabel('Intensity (a.u.)');
xlim([-2, 25]); ylim([0, 0.35]);

% Storage for report (Cell array to handle variable number of peaks)
% We will store matrices inside cells: [q_pos, Intensity, d_spacing]
peaks_vert_all = cell(1,6);
peaks_horz_all = cell(1,6);

% Define search range (Wider now to catch 3 peaks)
q_min = 2; 
q_max = 25;

% --- PEAK FINDING SETTINGS (TUNE THESE IF NEEDED) ---
min_prominence = 0.0001; % How much "sticking out" a peak must have (filters noise)
min_distance = 5.0;     % Min distance between peaks (filters dip-splitting)

for i = 1:6
    Current_Image = Im_corr{i};
    legend_list{i} = ['Step ' num2str(i)];
    
    % --- 1. EXTRACT DATA ---
    % Horizontal
    raw_h = interp2(XI, YI, Current_Image, x_line_horz, y_line_horz);
    % Vertical
    raw_v = interp2(XI, YI, Current_Image, x_line_vert, y_line_vert);
    
    % --- 2. GAP FILLING & SMOOTHING (CRITICAL FOR DIPS) ---
    % Step A: Moving Median (removes the sharp drop of the gap)
    gap_filled_h = movmedian(raw_h, 15, 'omitnan'); 
    gap_filled_v = movmedian(raw_v, 15, 'omitnan');
    
    % Step B: Gaussian Smooth (makes the hill smooth for peak finder)
    smooth_h = smoothdata(gap_filled_h, 'gaussian', 20);
    smooth_v = smoothdata(gap_filled_v, 'gaussian', 20);
    
    % --- 3. FIND PEAKS (Vertical) ---
    % Restrict to relevant q-range
    mask_v = q_axis_vert > q_min & q_axis_vert < q_max;
    q_roi_v = q_axis_vert(mask_v);
    I_roi_v = smooth_v(mask_v);
    
    [pks_v, locs_v] = findpeaks(I_roi_v, q_roi_v, ...
        'NPeaks', 3, ...              % Look for top 3 peaks
        'SortStr', 'descend', ...     % sort by height
        'MinPeakProminence', min_prominence, ... 
        'MinPeakDistance', min_distance); 
        
    % Store data: [q, I, d]
    if ~isempty(pks_v)
        peaks_vert_all{i} = [locs_v, pks_v, (2*pi)./locs_v];
    end

    % --- 4. FIND PEAKS (Horizontal) ---
    mask_h = q_axis_horz > q_min & q_axis_horz < q_max;
    q_roi_h = q_axis_horz(mask_h);
    I_roi_h = smooth_h(mask_h);
    
    [pks_h, locs_h] = findpeaks(I_roi_h, q_roi_h, ...
        'NPeaks', 3, ...
        'SortStr', 'descend', ...
        'MinPeakProminence', min_prominence, ...
        'MinPeakDistance', min_distance);
        
    if ~isempty(pks_h)
        peaks_horz_all{i} = [locs_h, pks_h, (2*pi)./locs_h];
    end
    
    % --- 5. PLOTTING ---
    % Plot the SMOOTHED curve (so user sees what the algorithm saw)
    plot(ax1, q_axis_vert, smooth_v, 'Color', colors(i,:), 'LineWidth', 1.5);
    plot(ax2, q_axis_horz, smooth_h, 'Color', colors(i,:), 'LineWidth', 1.5);
    
    % Mark the peaks
    if ~isempty(pks_v)
        plot(ax1, locs_v, pks_v, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
        % Label the secondary peaks (optional, uncomment if messy)
        % text(ax1, locs_v(1), pks_v(1)*1.05, '\downarrow', 'Color', colors(i,:), 'HorizontalAlignment', 'center');
    end
    if ~isempty(pks_h)
        plot(ax2, locs_h, pks_h, 'o', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
end

% Add Legends
legend(ax1, legend_list, 'Location', 'northeast');
legend(ax2, legend_list, 'Location', 'northeast');

% Save Figure
exportgraphics(figure(3), 'WAXS_MultiPeaks.png', 'Resolution', 300);

%% --- 7. AZIMUTHAL POLAR PLOTS (Circular Intensity Profile) ---
% This creates a "radar" or "polar" style plot.
% The distance from the center represents Intensity.
% Question 7: Look at variation around a circle at constant q.
% We need to find the q-radius of the main peak first.

% 2. Define Parameters
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
    peak_anisotropy(i) = max(smooth_I) - min(smooth_I);
end

legend('Location', 'northeastoutside'); % Move legend out of the way

% Optional: Adjust limits to make the shape clearer
% If the center is empty, you can set rlim([0 max_intensity])
rlim([0, inf]);

% Save Figure
exportgraphics(figure(4), 'WAXS_azimutha6.png', 'Resolution', 300);

%% --- 7. AZIMUTHAL POLAR PLOTS (Circular Intensity Profile) ---
% This creates a "radar" or "polar" style plot.
% The distance from the center represents Intensity.
% Question 7: Look at variation around a circle at constant q.
% We need to find the q-radius of the main peak first.

% 2. Define Parameters
target_q = 12.08;  % PEAK POSITION 
q_tol = 0.008;       % Width of the integration ring (0.04 for Peak1, 0.4 Peak 2)
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
    peak_anisotropy(i) = max(smooth_I) - min(smooth_I);
end

legend('Location', 'northeastoutside'); % Move legend out of the way

% Optional: Adjust limits to make the shape clearer
% If the center is empty, you can set rlim([0 max_intensity])
rlim([0, inf]);

% Save Figure
exportgraphics(figure(4), 'WAXS_azimutha12.png', 'Resolution', 300);
%% --- RELATE TO MACROSCOPIC CURVE ---
% Question 8: Relate observations to the Force-Displacement curve.

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
% We assume xpoints has 6 rows corresponding to the 6 images analyzed.
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

%% --- 8. PEAK AMPLITUDE VS AZIMUTH (Cartesian Plot) ---
% Objective: Quantify anisotropy by plotting Peak Intensity vs Angle.

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
            % Note: If data is very noisy, using 'quantile(vals, 0.98)' is safer than 'max'
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

%% --- 9. PEAK POSITION (Q-VALUE) VS AZIMUTH ---
% Objective: See if the inter-chain spacing (d-spacing) changes with angle.

figure(6); clf;
set(gcf, 'Position', [100, 100, 1000, 600]); 
hold on; box on;
title('Structural Distortion: Peak Position (q) vs. Azimuth');
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
            % Note: "max" is simple but can be noisy (pixel quantization). 
            % If the curve is too jumpy, we will smooth it later.
            [~, max_idx] = max(I_vals);
            peak_q_vals(k) = q_vals(max_idx);
        end
    end
    
    % 4. Smooth the curve
    % q-values are "quantized" by pixel size, so the raw data looks like a staircase.
    % Heavy smoothing is needed to see the real physical trend.
    peak_q_smooth = smoothdata(peak_q_vals, 'gaussian', 10);
    
    % 5. Plot
    plot(sector_deg, peak_q_smooth, 'Color', colors(i,:), ...
        'LineWidth', 2, 'DisplayName', ['Step ' num2str(i)]);
end

legend('Location', 'northeastoutside');
exportgraphics(figure(6), 'WAXS_Peak_Position_vs_Azimuth.png', 'Resolution', 300);

%% --- 9. PEAK POSITION (Q-VALUE) VS AZIMUTH ---
% Objective: See if the inter-chain spacing (d-spacing) changes with angle.

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
            % Note: "max" is simple but can be noisy (pixel quantization). 
            % If the curve is too jumpy, we will smooth it later.
            [~, max_idx] = max(I_vals);
            peak_q_vals(k) = q_vals(max_idx);
        end
    end
    
    % 4. Smooth the curve
    % q-values are "quantized" by pixel size, so the raw data looks like a staircase.
    % Heavy smoothing is needed to see the real physical trend.
    peak_q_smooth = smoothdata(peak_q_vals, 'gaussian', 10);
    
    % 5. Plot
    plot(sector_deg, peak_q_smooth, 'Color', colors(i,:), ...
        'LineWidth', 2, 'DisplayName', ['Step ' num2str(i)]);
end

legend('Location', 'northeastoutside');
exportgraphics(figure(7), 'WAXS_Peak_Position_vs_Azimuth.png', 'Resolution', 300);

