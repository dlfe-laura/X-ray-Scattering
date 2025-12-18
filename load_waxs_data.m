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
end

figure(2)
for i = 1:6
    ax = subplot(2, 3, i);
    pcolor(ax, Im_corr{i})      % Plotting the corrected data
    shading interp
    daspect([1 1 1])
    title(['Corrected Step ', num2str(i)])
    colorbar                    % Adds a color scale bar to see intensity values
end

%% EXTRACT 1D LINE PROFILES

% Get the sizes of the images
I1 = Im{1}; 
[sizex, sizey] = size(I1);

% Create a matrix of x and y coordinates
[XI, YI] = meshgrid(1:sizey,1:sizex);

% Determine the points defining the line
line = [0:sizex-centre_x];
x = centre_x + line;
y = ones(1,length(x))*centre_y;

% Interpolate the data from the 2D grid
I_line = interp2(XI,YI,I1,x,y);
