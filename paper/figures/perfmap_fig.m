% Postprocess the performance map figure created in the synthesis toolbox
% 
% This script is partly based on case_study/ParRob_nullspace_trajectory.m 
% from https://github.com/SchapplM/robotics-paper_icinco2021
% (ICINCO paper "Singularity Avoidance of Task-Redundant Robots in Pointing
% Tasks [...]", Schappler et al. 2021)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-01
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear

%% Initialization
importdir = ark2022_3T1R_dimsynth_data_dir();
datadir = fullfile(fileparts(which('ark2022_3T1R_dimsynth_data_dir.m')),'data');
tmp = load(fullfile(datadir, 'robot_groups.mat'));
paperfigdir = fullfile(datadir, '..', 'paper', 'figures');
RobotGroups = tmp.RobotGroups;

%% Create one performance map for each robot in the results set
for i = [1 7]%7%[]% 1:size(RobotGroups,1)
  GroupName = RobotGroups{i,1};
  if RobotGroups{i,3} == 0, continue; end % keine Ergebnisse vorliegend
  fprintf('Zeichne Redundanzkarte für PKM-Gruppe %d/%d (%s)\n', i, size(RobotGroups,1), GroupName);
  erg = load(fullfile(datadir, sprintf('detail_result_group_%s.mat', GroupName)));
  tmpdir_i = fullfile(importdir, erg.OptName, 'tmp', sprintf('%d_%s', ...
    erg.LfdNr, erg.RobName));
  fprintf('Ergebnis-Ordner: %s\n', tmpdir_i);
  % Einstellungen generieren
  setfile = dir(fullfile(importdir, erg.OptName, '*settings.mat'));
  d1 = load(fullfile(importdir, erg.OptName, setfile(1).name));
  Set_i = cds_settings_update(d1.Set);
  [R, Structure_i] = cds_dimsynth_robot(Set_i, d1.Traj, d1.Structures{erg.LfdNr}, true);

  trajdatafiles = dir(fullfile(tmpdir_i, '*_Traj*.mat'));
  perfmapfiles = dir(fullfile(tmpdir_i, '*Konfig*TaskRedPerfMap_Data.mat'));
  trajstats = array2table(NaN(length(trajdatafiles),4), 'VariableNames', ...
    {'ConfigNum', 'perfmapfileidx', 'Fval', 'TrajNum'});

  wn_sel = zeros(R.idx_ik_length.wntraj,1);
  PHIz_traj = NaN(length(d1.Traj.t), length(trajdatafiles));
  for k = 1:length(trajdatafiles)
    tmp = load(fullfile(tmpdir_i, trajdatafiles(k).name));
    [tokens,~] = regexp(trajdatafiles(k).name, 'Konfig(\d)+', 'tokens', 'match');
    confignum = str2double(tokens{1}{1});
    perfmapfileidx = find(contains({perfmapfiles.name}, sprintf('Konfig%d', confignum)));
    row_k = array2table([confignum, perfmapfileidx, tmp.fval, tmp.i_ar-1]);
    row_k.Properties.VariableNames = trajstats.Properties.VariableNames;
    trajstats(k,:) = row_k;
    wn_sel = wn_sel + tmp.s.wn;
    PHIz_traj(:,k) = tmp.X2(:,6);
  end
  k_iO = trajstats.Fval <= 1e3;
  k_plot = find(k_iO, 1, 'first');
  % Wähle weitere i.O.-Trajektorien dieser Konfiguration
  k_iOc = k_iO & trajstats.ConfigNum == trajstats.ConfigNum(k_plot);
  fprintf('%d/%d Trajektorien führen zu erfolgreichem Ergebnis. %d für diese Konfiguration.\n', ...
    sum(k_iO), length(trajdatafiles), sum(k_iOc));

  % Redundanzkarte laden (passend zur Nummer der Konfiguration)
  if isempty(perfmapfiles), error('Datei nicht gefunden'); end
  dpm = load(fullfile(tmpdir_i, perfmapfiles(trajstats.perfmapfileidx(k_plot)).name));

  wn_sel(:) = 0;
  wn_sel(R.idx_iktraj_wnP.jac_cond) = 1;
%   wn_sel(R.idx_iktraj_wnP.ikjac_cond) = 1;

  jj=0; I_wn_traj = zeros(R.idx_ik_length.wnpos,1);
  for f = fields(R.idx_ikpos_wn)'
    jj=jj+1; I_wn_traj(jj) = R.idx_iktraj_wnP.(f{1});
  end
  % Bild plotten (ohne Trajektorien)
  cds_debug_taskred_perfmap(Set_i, Structure_i, dpm.H_all, dpm.s_ref, ...
    dpm.s_tref, dpm.phiz_range, NaN(length(dpm.s_tref),1), NaN(length(dpm.s_tref),1), ...
    struct('wn', wn_sel(I_wn_traj), ...
    'i_ar', 0, 'name_prefix_ardbg', '', 'fval', 0, 'logscale', true, ...
    'critnames', {fields(R.idx_ikpos_wn)'}, 'constrvioltext', ''));

  % Zeichne Trajektorien ein (die für diese Konfiguration gelten)
  for k = find(k_iOc)'
    trajhdl = plot(dpm.s_tref, 180/pi*PHIz_traj(:,k), 'k-', 'LineWidth', 2);
    set(trajhdl, 'DisplayName', 'Traj')
  end
  sgtitle(''); title(''); % erst hiernach children bestimmen

  fighdl = gcf();
  fch = get(fighdl, 'children');
  axhdl = fch(strcmp(get(fch, 'type'), 'axes'));
  axch = get(axhdl, 'children');
  legdummyhdl = [];
  leglbl = {};
  for jj = 1:length(axch)
    set(axch(jj), 'MarkerSize', 1)
    if strcmp(get(axch(jj), 'DisplayName'), 'Joint Lim')
      delete(axch(jj)); % Gelenkgrenzen hier nicht verwertbar
      continue
    end
    if strcmp(get(axch(jj), 'DisplayName'), 'Act. Sing.')
      set(axch(jj), 'DisplayName', 'Sing. Type II')
    end
    if strcmp(get(axch(jj), 'DisplayName'), 'IK Sing.')
      set(axch(jj), 'DisplayName', 'Singularity Type I')
    end
    if ~strcmp(get(axch(jj), 'Marker'), 'none')
      leglbl = [leglbl, get(axch(jj), 'DisplayName')];
      legdummyhdl = [legdummyhdl; plot(NaN,NaN)]; %#ok<AGROW> 
      set(legdummyhdl(end), 'Marker', get(axch(jj), 'Marker'));
      set(legdummyhdl(end), 'MarkerSize', 5);
      set(legdummyhdl(end), 'DisplayName', get(axch(jj), 'DisplayName'));
      set(legdummyhdl(end), 'Color', get(axch(jj), 'Color'));
      set(legdummyhdl(end), 'LineStyle', get(axch(jj), 'LineStyle'));
    end
  end
  legdummyhdl(strcmp(get(legdummyhdl, 'DisplayName'), 'Traj')) = trajhdl;
  
  % Bestimme Grenzen des Datenbereichs
%   Hcond = dpm.H_all(:,:,R.idx_ikpos_wn.jac_cond);
%   all(isnan(Hcond),2)
  ylim(180/pi*minmax2(PHIz_traj(:,k_plot)')+[-45, +45])

  % Speichere das Bild im Paper-Format
  figure_format_publication(fighdl);
  set(gca, 'xticklabel', {});
  set_size_plot_subplot(fighdl, ...
    11.6,4,axhdl,...
    0.09,0.16,0.14,0.12,... %l r u d
    0,0) % x y
  xlabel('Normalized trajectory progress $s$', 'interpreter', 'latex');
  ylabel('Red. coord. $\varphi_z$ in deg', 'interpreter', 'latex');
  cbhdl = fch(strcmp(get(fch, 'type'), 'colorbar'));
  set(cbhdl, 'Ticks', [1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6], 'TickLabels', ...
    {'1', '10', '10^{2}', '10^{3}', '10^{4}', '10^{5}', '10^{6}'})
  cbyh = ylabel(cbhdl,'Performance criterion (cond.)', ...
    'Rotation',90, 'interpreter', 'latex');
%   [leghdl, objsb] = legend(legdummyhdl);
%   leghdl = fch(strcmp(get(fch, 'type'), 'legend'));
%   set(leghdl, 'orientation', 'horizontal', 'position', [0.15,0.92,0.7,0.05]);
%   delete(leghdl)
  leghdl = legendflex(legdummyhdl, leglbl, 'anchor', {'n','n'}, ...
    'ref', fighdl, ... % an Figure ausrichten (mitten oben)
    'buffer', [0 -1], ... % Kein Versatz notwendig, da mittig oben
    'ncol', 0, 'nrow', 1, ... % eine Zeile für Legende
    'fontsize', 8, ...
    'xscale', 0.6, ... % Kleine Symbole
    'padding', [-3,-3,3], ... % Leerraum reduzieren
    'box', 'on');
  % Schriftart der Legende aktualisieren
  figure_format_publication(fighdl);
  drawnow();

  t1 = tic();
  % For this to work, the Java heap memory has to be high enough.
  exportgraphics(fighdl, fullfile(paperfigdir, sprintf( ...
    'group%d_perfmap.pdf', i)),'ContentType','vector');
  saveas(fighdl, fullfile(paperfigdir, sprintf( ...
    'group%d_perfmap.fig', i)));
  fprintf('Exported performance map as vector graphics. Duration: %1.1fs\n', toc(t1));
  return
  % Create smaller version of the figures for manually stitching the
  % overall image
  xlim_orig = xlim();
  xlim(xlim_orig/2);
  set_size_plot_subplot(fighdl, ...
    4.7,3,axhdl,... % measure width in InkScape to obtain the value
    0.12,0.02,0.14,0.12,... %l r u d
    0,0) % x y
  delete(get(axhdl, 'xlabel'));
  delete(get(axhdl, 'ylabel'));
  delete(leghdl);
  delete(cbhdl);
  drawnow();
  exportgraphics(fighdl, fullfile(paperfigdir, sprintf( ...
    'group%d_perfmap_cut.pdf', i)),'ContentType','vector');
  saveas(fighdl, fullfile(paperfigdir, sprintf( ...
    'group%d_perfmap_cut.fig', i)));
  fprintf('Exported performance map without addition texts as vector graphics. Duration: %1.1fs\n', toc(t1));
end
return
%% Create one performance map figure for the document
% Did not work.
idx_sel = [1 7];
filenames = cell(1,2);
for i = 1:2
  filenames{i} = fullfile(fullfile(paperfigdir, sprintf( ...
      'group%d_perfmap.fig', idx_sel(i))));
end
% Erstes Bild öffnen und zweites hineinkopieren
close all
uiopen(filenames{1}, 1)
fig1hdl = gcf;
uiopen(filenames{2}, 1)
fig2hdl = gcf;
% Kopiere ax-Handle von Fig 2 nach Fig 1
fig1ch = get(fig1hdl,'children');
fig2ch = get(fig2hdl,'children');
ax1hdl = fig1ch(strcmp(get(fig1ch,'type'),'axes'));
ax2hdl = fig1ch(strcmp(get(fig2ch,'type'),'axes'));
% delete(ax1hdl(2))
copyobj(ax2hdl(2), fig1hdl);
close(fig2hdl);
fig1ch_neu = get(fig1hdl,'children');
ax1hdl_neu = fig1ch_neu(strcmp(get(fig1ch_neu,'type'),'axes'));
% Plot neu formatieren
% delete(ax1hdl_neu(1));
get(ax1hdl_neu(2), 'position')
get(ax1hdl_neu(3), 'position')
set_size_plot_subplot(fig1hdl, ...
  11.6,4,ax1hdl_neu(2:3)',...
  0.09,0.16,0.14,0.12,... %l r u d
  0.05,0) % x y
drawnow()
