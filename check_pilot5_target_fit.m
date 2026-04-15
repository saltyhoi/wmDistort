%% check_pilot5_target_fit.m
% Check whether all target positions fall inside each pilot5 aperture

clear; clc;

%% -----------------------------
% User-set parameters
%% -----------------------------
screen_height_px    = 1440;   % example
screen_height_cm    = 30;     % behavioral room
viewing_distance_cm = 52;     % behavioral room

wm_ecc      = 6.5;              % deg
aperture_size = 15;           % kept here for reference
phi         = 0.618;

repetitions = 16;
start_deg   = 11.25;
x_offset_deg = 4;

diag_scale = 0.85;            % same scale you plan to use for diagonals

%% -----------------------------
% Derived quantities
%% -----------------------------
ppd = screen_height_px / (2 * atan2d(screen_height_cm/2, viewing_distance_cm));
center = [0, 0];  % use aperture-centered coordinates for the math

degrees = start_deg + (0:(repetitions-1)) * (360/repetitions);

% IMPORTANT:
% In your codebase, these names are visually swapped by convention:
% aperture_wide_rect_px = [0 0 H*phi H];  % actually tall-looking
% aperture_tall_rect_px = [0 0 H H*phi];  % actually wide-looking
H = screen_height_px;

% For left/right wide rectangles, use the same dimensions as your task code
rect_w_lr = H ;
rect_h_lr = H * phi;

% For diagonal rectangles, use the scaled dimensions
rect_w_diag = H * diag_scale;
rect_h_diag = H * phi * diag_scale;

x_offset_px = x_offset_deg * ppd;

%% -----------------------------
% Generate target positions
%% -----------------------------
target_xy_deg = wm_ecc * [cosd(degrees(:)), sind(degrees(:))];
target_xy_px  = ppd * target_xy_deg;   % [x y] in centered coords

%% -----------------------------
% Check each aperture
%% -----------------------------
results = struct();

% left/right shifted axis-aligned rectangles
results.left_wide_rect  = checkAxisAlignedRect(target_xy_px, center + [-x_offset_px, 0], rect_w_lr, rect_h_lr);
results.right_wide_rect = checkAxisAlignedRect(target_xy_px, center + [ x_offset_px, 0], rect_w_lr, rect_h_lr);

% rotated rectangles
results.diag45_rect  = checkRotRect(target_xy_px, center, rect_w_diag, rect_h_diag, 45);
results.diag135_rect = checkRotRect(target_xy_px, center, rect_w_diag, rect_h_diag, 135);

%% -----------------------------
% Print summary
%% -----------------------------
shape_names = fieldnames(results);

fprintf('\n=== Pilot5 target-fit check ===\n');
fprintf('ppd = %.3f px/deg\n', ppd);
fprintf('wm_ecc = %.3f deg (%.3f px)\n', wm_ecc, wm_ecc * ppd);
fprintf('x_offset_deg = %.3f (%.3f px)\n', x_offset_deg, x_offset_px);
fprintf('diag_scale = %.3f\n\n', diag_scale);

for i = 1:numel(shape_names)
    shape = shape_names{i};
    res = results.(shape);

    fprintf('%s\n', shape);
    fprintf('  inside count: %d / %d\n', sum(res.inside), numel(res.inside));

    if all(res.inside)
        fprintf('  status: PASS\n\n');
    else
        fprintf('  status: FAIL\n');
        fprintf('  failed angles (deg): ');
        fprintf('%.2f ', degrees(~res.inside));
        fprintf('\n\n');
    end
end

%% -----------------------------
% Optional plot
%% -----------------------------
figure; hold on; axis equal;
title('Pilot5 aperture check');
xlabel('x (px)');
ylabel('y (px)');

% plot target points
plot(target_xy_px(:,1), target_xy_px(:,2), 'ko', 'MarkerFaceColor', 'k');

% draw apertures
drawAxisAlignedRect(center + [-x_offset_px, 0], rect_w_lr, rect_h_lr, 'b');
drawAxisAlignedRect(center + [ x_offset_px, 0], rect_w_lr, rect_h_lr, 'r');
drawRotRect(center, rect_w_diag, rect_h_diag, 45,  'm');
drawRotRect(center, rect_w_diag, rect_h_diag, 135, 'g');

legend({'Targets','Left wide','Right wide','Diag45','Diag135'});
hold off;

%% =============================
% Local functions
%% =============================

function out = checkAxisAlignedRect(points_xy, rect_center, w, h)
    x = points_xy(:,1) - rect_center(1);
    y = points_xy(:,2) - rect_center(2);

    inside = (abs(x) <= w/2) & (abs(y) <= h/2);

    out.inside = inside;
end

function out = checkRotRect(points_xy, rect_center, w, h, angle_deg)
    theta = deg2rad(angle_deg);

    % rotate points into rectangle coordinates
    R = [cos(theta)  sin(theta);
        -sin(theta)  cos(theta)];

    shifted = points_xy - rect_center;
    rotated = (R * shifted')';

    inside = (abs(rotated(:,1)) <= w/2) & (abs(rotated(:,2)) <= h/2);

    out.inside = inside;
end

function drawAxisAlignedRect(rect_center, w, h, colorSpec)
    x0 = rect_center(1);
    y0 = rect_center(2);

    corners = [x0-w/2, y0-h/2;
               x0+w/2, y0-h/2;
               x0+w/2, y0+h/2;
               x0-w/2, y0+h/2;
               x0-w/2, y0-h/2];
    plot(corners(:,1), corners(:,2), colorSpec, 'LineWidth', 1.5);
end

function drawRotRect(rect_center, w, h, angle_deg, colorSpec)
    x0 = rect_center(1);
    y0 = rect_center(2);

    base = [-w/2, -h/2;
             w/2, -h/2;
             w/2,  h/2;
            -w/2,  h/2;
            -w/2, -h/2]';

    theta = deg2rad(angle_deg);
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

    rot = R * base;
    rot(1,:) = rot(1,:) + x0;
    rot(2,:) = rot(2,:) + y0;

    plot(rot(1,:), rot(2,:), colorSpec, 'LineWidth', 1.5);
end