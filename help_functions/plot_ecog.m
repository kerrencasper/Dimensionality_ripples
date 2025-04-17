function plot_ecog(plot_vec, mesh_pth, epos, limits, transparency, view_pos, elec_size)
% Modified plot_ecog function to plot electrodes in MNI space with specific colors
% assigned to each contact.
% @plot_vec now contains an nx3 matrix of RGB values for each contact.

% Load the mesh files
ld = load(fullfile(mesh_pth, 'surface_pial_right.mat'));
mesh_r = ld.mesh;

ld = load(fullfile(mesh_pth, 'surface_pial_left.mat'));
mesh_l = ld.mesh;

% Set default values if necessary
if nargin < 7
    elec_size = 20;
end
if nargin < 6
    view_pos = [-90 0];
end
if nargin < 5
    transparency = 0.5; % transparency of the mesh
end
if nargin < 4
    limits = [0 1]; % Default limits for consistency
end

% Plot the brain meshes with specified transparency
ft_plot_mesh(mesh_l, 'facealpha', transparency);
hold on;
ft_plot_mesh(mesh_r, 'facealpha', transparency);
lighting gouraud;
camlight;

% Plot the electrode positions with the specified RGB colors
ft_plot_mesh(epos, 'vertexcolor', plot_vec, 'vertexsize', elec_size);

% Set the color limits and colorbar
caxis(limits);
colorbar off; % Disable colorbar since RGB colors are used directly

% Set figure background to white
set(gcf, 'color', 'w');

% Set the specified view for the plot
view(view_pos);
end
