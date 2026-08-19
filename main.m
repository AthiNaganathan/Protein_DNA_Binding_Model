%% Initialization

% Check correct working directory
if ~isfile(fullfile(pwd, "binding_model.m"))
    error("Incorrect working directory, cannot find ./binding_model.m/")
end
mkdir(fullfile(pwd, "results")) % Create ./results/ folder, if missing

clear
clc

global T;
global s_min;
global concs;
global exp_S;
global err_S;

bm = binding_model;

%% Input parameters

pname = "NHP6A";        % Filename containing protein information in ./data/
n = 15;                 % No. of base pairs of DNA
T = 298;                % Temperature of experiment (K)
sites = 7;              % No. of consecutive DNA base pairs that protein binds to
modes = 1;              % No. of binding modes that the protein has (either 1 or 2)

%% Verify input parameters

if n >= 3*sites
    error("More than two proteins binding to a single fragment is currently not supported. Kindly set n < 3*sites .")
end

if modes ~= 1 && modes ~= 2
    error("Only 1 and 2 binding modes are currently supported. Kindly set modes = 1 or 2.")
end

% Protein data should be present in a .dat file with columns -
% 1. Concentration (micromolar)
% 2. Experimentally measured anisotropy
% 3. Standard deviation of measurements
% File to be placed in ./data/

if ~isfile(fullfile(pwd, "data", pname + ".dat"))   % Check if input data is present
    error("Input data file for given protein name not found in ./data/")
end

data = readmatrix(fullfile(pwd, "data", pname + ".dat"));
if shape(data, 2) ~= 3
    error("Input data does not contain exactly 3 columns. Kindly check if correct data is provided.")
end
concs = data(:, 1).';   % Experimental concentrations (u mol)
exp_S = data(:, 2).';   % Measured anisotropies
err_S = data(:, 3).';   % Standard deviations of experimental data

%% Initialize parameters for fitting algorithm

s_min = mean(exp_S(1:5));       % Anisotropy @ 0 concentration is fixed
if modes == 1               % Need to fit [binding energy (B. E.), anisotropy at saturation]
    lb = [-50, 0];              % Lower bounds
    ub = [0, 0.4];              % Upper bounds 
    ival = [-29, 0.212];        % Initial guess
else                       % Need to fit [major B.E., minor B.E., anisotropy at saturation]
    lb = [-50, -50, 0];         % Lower bounds
    ub = [0, 0, 0.4];           % Upper bounds
    ival = [-35, -29, 0.212];   % Initial guess
end

%% Run optimization

fix_vals = [n modes sites];
[fval, ~, residual, ~, ~, ~, jacobian] = lsqcurvefit(@(vars, concs) bm.init_calc_sigs([vars(1:modes+1) fix_vals], concs, modes, sites), ival, concs, exp_S, lb, ub);
conf = nlparci(fval, residual, 'jacobian', jacobian);

perms = bm.generate_perms(1:modes, floor(n/sites));
if modes == 2
    perms(perms == 11) = [];    % Remove double strong binding
end
DNA_states = bm.fill_perms(perms, n, modes, sites);
[no_p, calc_S] = bm.get_calc_sigs(fval, concs, DNA_states, DNA_states, modes, sites);

% Get probabilities of only one or only two proteins binding
if modes == 1
    states_1 = bm.fill_perms(1, n, modes, sites); 
    [no_p_1, calc_S_1] = bm.get_calc_sigs(fval, concs, states_1, DNA_states, modes, sites);
    states_2 = bm.fill_perms(11, n, modes, sites); 
    [no_p_2, calc_S_2] = bm.get_calc_sigs(fval, concs, states_2, DNA_states, modes, sites);
else
    states_1 = bm.fill_perms([1, 2], n, modes, sites); 
    [no_p_1, calc_S_1] = bm.get_calc_sigs(fval, concs, states_1, DNA_states, modes, sites);
    states_2 = bm.fill_perms([12, 22], n, modes, sites); 
    [no_p_2, calc_S_2] = bm.get_calc_sigs(fval, concs, states_2, DNA_states, modes, sites);
end

%% Plot results
clf('reset')

c = newline;
text_b = "Name : " + pname + c + "No of sites : " + num2str(sites) + c + "Binding modes : " + num2str(modes);
text_n = "S min : " + s_min ;
if modes == 1
    text_n = text_n + c + "S max : " + num2str(round(fval(2),4)) + " +- " + num2str(round(max(abs(repelem(fval(2), 2)-conf(2, :))),4));
    text_n = text_n + c + "Binding f.e : " + num2str(round(fval(1),4)) + " +- " + num2str(round(max(abs(repelem(fval(1), 2)-conf(1, :))),4));
else
    text_n = text_n + c + "S max : " + num2str(round(fval(3),4)) + " +- " + num2str(round(max(abs(repelem(fval(3), 2)-conf(3, :))),4));
    text_n = text_n + c + "Major b.f.e : " + num2str(round(fval(1),4)) + " +- " + num2str(round(max(abs(repelem(fval(1), 2)-conf(1, :))),4));
    text_n = text_n + c + "Minor b.f.e : " + num2str(round(fval(2),4)) + " +- " + num2str(round(max(abs(repelem(fval(2), 2)-conf(2, :))),4));
end

t = tiledlayout(2, 2);

ax1 = nexttile([2 1]);
hold on;
h = errorbar(concs, exp_S, err_S, 'or');
set(get(h, 'Parent'), 'XScale', 'log');
semilogx(concs, calc_S, 'b', concs, calc_S_1, '--', concs, calc_S_2, '--');
xlabel('Concentration (u mol/L)')
ylabel('Anisotropy')
s_max = fval(end);
yline([s_min+(s_max-s_min) s_min+2*(s_max-s_min)], '--', {'1p saturating sig', '2p saturating sig'})

lgd = legend( "Experimental S", "Calculated S", "S from 1 prot(s) bound", "S from 2 prot(s) bound", '', '', 'location', 'northwest');
drawnow;
pos = lgd.Position;
lineHeight = 0.03;

ann_b = annotation('textbox', [pos(1), pos(2) - 1.5*lineHeight, pos(3), lineHeight], ...
    'String', text_b, 'FontWeight', 'bold', 'FontSize', 12, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top');
ann_n = annotation('textbox', [pos(1), pos(2) - 4*lineHeight, pos(3), lineHeight], ...
    'String', text_n, 'FontSize', 10, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top');
disp(text_b)
disp(text_n)

ax2 = nexttile;
fig2 = semilogx(concs, no_p);
title('Occupancy vs Conc')
xlabel('Concentration mol/L')
ylabel('Occupancy')

ax3 = nexttile;
[~, b_dis] = bm.get_b_dis(fval, concs, DNA_states, modes, sites);
h = surf(1:n, concs, b_dis);
set(get(h, 'Parent'), 'YScale', 'log')
title('Binding distributions')
xlabel('Base Pair')
ylabel('Concentration mol/L')
zlabel('Probability of binding')

exp_labels = {'Conc(u mol L-1)' 'Calculated S' 'Signal change for 1 protein' 'Signal change for 2 proteins' 'No. of proteins bound'};
exp_matrix = num2cell([concs' calc_S' calc_S_1' calc_S_2' no_p']);
exp = [exp_labels; exp_matrix];

%% Save calculated distributions

op_loc = fullfile(pwd, "results", string(pname));
mkdir(op_loc)
rmdir(op_loc, 's')
mkdir(op_loc)

% Save signal changes per protein bound as .xlsx
delete(fullfile(op_loc, string(pname)+".xlsx"))
writecell(exp, fullfile(op_loc, string(pname)+".xlsx"))

% Save the displayed graph
saveas(gcf, fullfile(op_loc, string(pname)+".png"))

% Save calculated anisotropy as .mat
p_var = {pname, concs, exp_S, err_S, calc_S};
save(fullfile(op_loc, string(pname) + ".mat"), "p_var");

fprintf("\nOutput files saved to %s\n", op_loc)

%% Done

disp("Done")