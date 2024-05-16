clear all
close all


           % Manually define data matrix data
data = dlmread("dem_output1.txt");

% Manually define data matrix data1
data1 = dlmread("Intrpolated_SMART_FARM1024x1024.txt")
                z = data;
                nbligne = numel(unique(data(:, 2))); % Number of unique y coordinates
                nbcolonne = numel(unique(data(:, 1))); % Number of unique x coordinates
                disp('Processed file with other than three columns');

zz=z;


   
   imin=1;
   jmin=1;
   imax=numel(unique(data1(:, 2))); 
   jmax=numel(unique(data1(:, 1)));
drained_areas = zeros(imax-imin+1, jmax-jmin+1);
%Loop over each pixel
for k = imin:imax
    for f = jmin:jmax
        % Determine the sub-watershed of the current pixel
      iy_exutoire = k;
      jx_exutoire = f;
    x_outlet = jx_exutoire;
y_outlet = iy_exutoire;
[M, N] = size(data1);
Nx =  M; % Or any calculation that reflects the dimensions of zz
Ny = N; % Same as above
A(1:Nx,1:Ny)=-1;                  % Indicator function A is initially a non-drained region...
A(jx_exutoire,iy_exutoire)=+1;        % ... except for a single pixel, the "outlet" pixel, 
Z_outlet=zz(jx_exutoire,iy_exutoire);  % ......which belongs initially to the drained region.
display(Z_outlet);

AA(1:Nx+2,1:Ny+2)=-1;           % Indicator function AA is initially a non-drained region
ZZ(1:Nx+2,1:Ny+2)=-Inf;         % Pixels outside domain have -Inf altitude
AA(2:Nx+1,2:Ny+1)=A(1:Nx,1:Ny); % Extended indicator function
ZZ(2:Nx+1,2:Ny+1)=zz(1:Nx,1:Ny); % Altitudes of all the pixels of AA

AA0=AA; % Initializing AA0 (RA24Mar2024)!

% ----------------------------------------------
% OptionSweepJI : sweep directions (RA24Mar2024)! 
% ----------------------------------------------
% OptionSweepJI=1 to sweep only as (J=...,I=...)
% OptionSweepJI=2 to sweep only as (I=...,J=...)
% OptionSweepJI=12 to sweep alternately as (J=...,I=...) then (I=...,J=...)
% -------------------------------------------------------------------------
  OptionSweepJI=2; 
% -------------------------------------------------------------------------

% -------------------------------------------------------
% Choose the maximum number of sweeps to do (RA24Mar2024)! 
% -------------------------------------------------------
OptionSweepJI=12; 
% -------------------------------------------------------------------------

% -------------------------------------------------------
% Choose the maximum number of sweeps to do (RA24Mar2024)! 
% -------------------------------------------------------
  Ksweeps = 20; % Example: Ksweeps=5 (RA24Mar2024)
% -------------------------------------------------------
ksweepCounter = 0; % Initialize the sweep counter
stability = true; % Assume stability initially
stabilityNumber=Inf; % Initialize this number to infinity
stabilityCounter = 0; % Counter for stable iterations

while ksweepCounter < Ksweeps && stabilityCounter <= 2
    % Increment the sweep counter:
      ksweepCounter = ksweepCounter + 1; % RA01Apr2024
    %
    % (( Copy the current state of AA to use in the next iteration ))
    % AA=AA0; % NO, NOT HERE (see end of sweep loop instead)! RA24Mar2024!
    %
    if OptionSweepJI==1 | OptionSweepJI==12, 
    for jy=1:Ny
        for ix=1:Nx
            iixx=1+ix;
            jjyy=1+jy;
             if AA(iixx,jjyy) == +1 % Test it only if A=+1
                if  ZZ(iixx+1,jjyy)   >= ZZ(iixx,jjyy), AA(iixx+1,jjyy)=+1; end
                if  ZZ(iixx-1,jjyy)   >= ZZ(iixx,jjyy), AA(iixx-1,jjyy)=+1; end
                if  ZZ(iixx,jjyy+1)   >= ZZ(iixx,jjyy), AA(iixx,jjyy+1)=+1; end
                if  ZZ(iixx,jjyy-1)   >= ZZ(iixx,jjyy), AA(iixx,jjyy-1)=+1; end
                if  ZZ(iixx+1,jjyy+1) >= ZZ(iixx,jjyy), AA(iixx+1,jjyy+1)=+1; end
                if  ZZ(iixx-1,jjyy+1) >= ZZ(iixx,jjyy), AA(iixx-1,jjyy+1)=+1; end
                if  ZZ(iixx+1,jjyy-1) >= ZZ(iixx,jjyy), AA(iixx+1,jjyy-1)=+1; end
                if  ZZ(iixx-1,jjyy-1) >= ZZ(iixx,jjyy), AA(iixx-1,jjyy-1)=+1; end
            end
        end
    end
    % AA0 = AA;   % Update AA0 here ? NON !
    end % End of first sweep direction
    %
    if OptionSweepJI==2 | OptionSweepJI==12, 
    for jy=1:Ny
        for ix=1:Nx
            iixx=1+ix;
            jjyy=1+jy;
            %% if AA0(iixx,jjyy) == +1 % Test it only if AA0=+1 <= NON ! RA24Mar2024
            if AA(iixx,jjyy) == +1 % Test it only if AA=+1
                if  ZZ(iixx+1,jjyy)   >= ZZ(iixx,jjyy), AA(iixx+1,jjyy)=+1; end
                if  ZZ(iixx-1,jjyy)   >= ZZ(iixx,jjyy), AA(iixx-1,jjyy)=+1; end
                if  ZZ(iixx,jjyy+1)   >= ZZ(iixx,jjyy), AA(iixx,jjyy+1)=+1; end
                if  ZZ(iixx,jjyy-1)   >= ZZ(iixx,jjyy), AA(iixx,jjyy-1)=+1; end
                if  ZZ(iixx+1,jjyy+1) >= ZZ(iixx,jjyy), AA(iixx+1,jjyy+1)=+1; end
                if  ZZ(iixx-1,jjyy+1) >= ZZ(iixx,jjyy), AA(iixx-1,jjyy+1)=+1; end
                if  ZZ(iixx+1,jjyy-1) >= ZZ(iixx,jjyy), AA(iixx+1,jjyy-1)=+1; end
                if  ZZ(iixx-1,jjyy-1) >= ZZ(iixx,jjyy), AA(iixx-1,jjyy-1)=+1; end
            
            end
        end
    end
    % AA0 = AA;   % Update AA0 here ? NON !
    end % End of 2nd sweep direction

% %   % Increment the sweep count % MOVED AT THE BEGINNING OF WHILE LOOP !
% %     ksweepCounter = ksweepCounter + 1;

    % Check for stability
    if isequal(AA, AA0)
        fprintf('Sweep %d: Algorithm is stable\n', ksweepCounter);
        stabilityCounter = stabilityCounter + 1;
        stabilityNumber = min(stabilityNumber, ksweepCounter); % Update stability number
        if ksweepCounter == stabilityNumber
            disp('----------------------');
            FirstStableSweep = num2str(ksweepCounter);
            FirstStableSweep = ['First stable sweep: ' FirstStableSweep];
            disp(FirstStableSweep);
            disp('----------------------');
        end
    else
        fprintf('Sweep %d: Algorithm is not stable\n', ksweepCounter);
        stabilityCounter = 0; % Reset stability counter if not stable
    end
 A(1:Nx,1:Ny) = AA(2:Nx+1,2:Ny+1);
 if ksweepCounter > Ksweeps   % if ksweepCounter >= Ksweeps
         break;
 end
    % 
    AA0=AA;
end
zzz=zz;
outlet_x = jx_exutoire;
outlet_y = iy_exutoire;
outlet_z = zzz(iy_exutoire, jx_exutoire); % Use the elevation value from the DEM
outlet_z1 = zzz(jx_exutoire, iy_exutoire);% for stem3
figure(121);
surf(zzz);
hold on;

% Plot DEM surface for drained area with a different color
% Define zzz_drained with the same dimensions as zzz
zzz_drained = NaN(size(zzz));

zzz_drained(A ~= 1) = NaN; % Mask non-drained area with NaN
 %zzz_drained(D ~= 1) = NaN; % Mask non-drained area with NaN
surf(zzz_drained, 'FaceColor', 'red', 'EdgeColor', 'none');

% Define the maximum elevation value (zmax) of the DEM
zmax = max(zzz(:))+0.5;

% Define the endpoint of the stick at the outlet coordinates (jx_exutoire, iy_exutoire) with elevation equal to zmax
stick_start = [iy_exutoire, jx_exutoire, zzz(jx_exutoire, iy_exutoire)];
stick_end = [iy_exutoire, jx_exutoire, zmax];

% Plot the stick representing the outlet
quiver3(stick_start(1), stick_start(2), stick_start(3), ...
    stick_end(1)-stick_start(1), stick_end(2)-stick_start(2), stick_end(3)-stick_start(3), ...
    'Color', 'yellow', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'Autoscale', 'off', 'ShowArrowHead', 'off', 'LineStyle', '-');

outlet_color = [0 0 0]; % Black color for the outlet pixel
scatter3(outlet_y, outlet_x, outlet_z1, 100, outlet_color, 'filled');

% Set axis labels, title, and other plot properties
xlabel('X Coordinate');
ylabel('Y Coordinate');
zlabel('Elevation');
title('DEM with Drained Area Highlighted');

deltax=1; 
deltay=1;
 K = (A+1) ./ 2;
        drained_area = sum(K(:)) * deltax * deltay; % Sum of all elements in K multiplied by the area of each pixel
      
        % Save the drained area corresponding to the current pixel
        drained_areas(k-imin+1, f-jmin+1) = drained_area;
        x_coords(k-imin+1, f-jmin+1) = f;
        y_coords(k-imin+1, f-jmin+1) = k;
   
% Save the drained areas matrix to a file
save('drained_areas.mat', 'drained_areas');
    end
end 
% Find the pixel that gives the maximum sub-watershed area
%}
[max_drained_area, idx] = max(drained_areas(:));
[k_max, f_max] = ind2sub(size(drained_areas), idx);
fprintf('Coordinates of pixel with the maximum drained area: (%d, %d)\n', k_max, f_max);
% Plot a 3D bar graph of the drained areas
figure;
bar3(drained_areas);
xlabel('X Coordinate');
ylabel('Y Coordinate');
zlabel('Drained Area');
title('3D Bar Graph of Drained Areas');

% Additionally, you may want to highlight the pixel with the maximum drained area
hold on;
plot3(f_max, k_max, max_drained_area, 'ro', 'MarkerSize', 10); % Plotting the pixel with the maximum drained area on the graph

% Add text annotations to each bar indicating the drained area and real coordinates
for i = 1:size(drained_areas, 1)
    for j = 1:size(drained_areas, 2)
        text(j, i, drained_areas(i, j), ['(', num2str(i), ',', num2str(j), ') - ', num2str(drained_areas(i, j))], 'FontSize', 8);
    end
end
hold off;
% Determine the coordinates of the outlet pixel on the surface of the DEM

