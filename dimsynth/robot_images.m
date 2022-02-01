% Create high-resolution images for the robots from the optimization
% 
% Preliminaries:
% * Aggregate data with eval_figures_pareto_groups.m
% * Assemble details with select_eval_robot_examples.m
% 
% This script is based on the same file eval_figures_pareto.m from 
% https://github.com/SchapplM/robsynth-paper_mhi2021
% (MHI paper "Combined Structural and Dimensional Synthesis of a Parallel
% Robot for Cryogenic Handling Tasks", Schappler et al. 2022)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-01
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear
close all

warning('off', 'Coder:MATLAB:rankDeficientMatrix'); % für invkin2
datadir = fullfile(fileparts(which('ark2022_3T1R_dimsynth_data_dir.m')),'data');
tmp = load(fullfile(datadir, 'robot_groups.mat'));
paperfigdir = fullfile(datadir, '..', 'paper', 'figures');
RobotGroups = tmp.RobotGroups;
namestablepath = fullfile(datadir, 'robot_names_latex.csv');
ResTab_NameTrans = readtable(namestablepath, 'Delimiter', ';');

%% Alle Gruppen durchgehen
countrob = 0; % Nummern aus Legende
for i = 1:size(RobotGroups,1)
  GroupName = RobotGroups{i,1};
  if RobotGroups{i,3} == 0, continue; end % keine Ergebnisse vorliegend
  fprintf('Zeichne Bild für PKM-Gruppe %d/%d (%s)\n', i, size(RobotGroups,1), GroupName);
  resfile = fullfile(datadir, sprintf('detail_result_group_%s.mat', GroupName));
  if ~exist(resfile, 'file')
    warning('Result file for Group %d (%s) does not exist: %s', i, GroupName, resfile);
    continue;
  end
  erg = load(resfile);
  % countrob = countrob + 1; % nur gültige PKM-Gruppen zählen
  countrob = i; % Geht nur, wenn alle Gruppen gültig sind.
  % Daten der Ergebnisse laden
  fprintf('Anzahl Variablen: %d\n', length(erg.Structure.vartypes));
  disp(strjoin(erg.Structure.varnames));
  fprintf('Parametertypen: [%s]\n', disp_array(erg.Structure.vartypes', '%1.0f'));
  close all
  %% Roboter initialisieren
  R = erg.R;
  % Dicke der Segmente reduzieren
  for kk = 1:R.NLEG
    R.Leg(kk).DesPar.seg_par(:,1) = 1e-3;
    R.Leg(kk).DesPar.seg_par(:,2) = 1e-2;
  end
  parroblib_addtopath({R.mdlname});
  %% Anpassung der Visualisierung
  Q = erg.Q;
  X = erg.X;
  Q_plot = Q;
  X_plot = X;

  %% Bild zeichnen und Formatieren
  change_current_figure(1); clf; hold on;
  s_plot = struct('ks_legs', [], 'ks_platform', [], 'mode', 4);
  % Plotte in mittlerer Stellung des Roboters
  x_mean = mean(minmax2(X_plot(:,1:3)'),2);
  % Finde den Punkt der Trajektorie, der der Mitte am nächsten ist.
  [~,I_mean] = min(sum((repmat(x_mean',size(X_plot,1),1)-X_plot(:,1:3)).^2,2));
  R.plot(Q_plot(I_mean,:)', X_plot(I_mean,:)', s_plot);
  view(3);
  title('');xlabel('');ylabel('');zlabel('');
  LegColors = [ [0 1 0]; [0 0 0]; [1 0 1]; [0 0 1]; [1 0 0] ]; % Grün, Schwarz, Violett, blau, rot
  % automatisch generierten Plot nachverarbeiten
  ch = get(gca, 'Children');
  for jj = 1:length(ch)
    % KS entfernen
    if strcmp(ch(jj).Type, 'hgtransform')
      delete(ch(jj)); continue
    end
%     % Dreieck entfernen, das die Basis anzeigt. Kann man besser in InkScape
%     % neu zeichnen
%     if strcmp(ch(jj).Type, 'line')
%       delete(ch(jj)); continue
%     end
    % Weise den Beinketten neue Farben zu
    for kk = 1:R.NLEG
      if contains(get(ch(jj), 'DisplayName'), sprintf('Leg_%d_Link', kk))
        set(ch(jj), 'EdgeColor', LegColors(kk,:));
      end
    end 
  end

  % Zeichne die Trajektorie (nur in einen Roboter)
  if i == 2
    Xt_W = NaN(size(X_plot,1),3);
    for k = 1:size(X_plot,1)
      Xt_W(k,:) = eye(3,4)*R.T_W_0 * [X_plot(k,1:3)';1];
    end
    plot3(Xt_W(:,1), Xt_W(:,2), Xt_W(:,3), 'c-', 'LineWidth', 2);
  end
  
  % Dummy-Plot für Legende mit Roboter-Marker
  hdl=plot3(0,0,NaN, RobotGroups{i,4});
  % abgespeicherten Marker-Typ einzeichnen (zur Wiedererkennung beim
  % Zusammenstellen der Bilder)
  % lh = legend(hdl, legstr_i, 'interpreter', 'latex'); 
  % set(lh, 'location', 'northeastoutside');
  
  % Passe die Perspektive an (wenig Verdeckungen)
  view(46,27); % Standard-Perspektive
  % Manuelle Anpassungen für einzelne Bilder möglichst wenig Verdeckungen
  if strcmp(GroupName, 'P4PRRRR4V1G') || strcmp(GroupName, 'P4RPRRR8V1G')
    view(45,16);
  elseif strcmp(GroupName, 'P4RRRRR7G') 
    view(100,30);
  end
  set(gca,'XTICKLABEL',{});set(gca,'YTICKLABEL', {});set(gca,'ZTICKLABEL',{});
  set(gca,'xtick',[],'ytick',[],'ztick',[]);
  set(get(gca, 'XAxis'), 'visible', 'off');
  set(get(gca, 'YAxis'), 'visible', 'off');
  set(get(gca, 'ZAxis'), 'visible', 'off');
  figure_format_publication(gca);
  set(gca, 'Box', 'off');
  set(1, 'windowstyle', 'normal');
  set_size_plot_subplot(1, ...
    8,8,gca,...
    0,0,0,0,0,0)
  drawnow();
  % Bild speichern
  name = sprintf('RobotFig_PlotNum%d_RobGroup%d_%s', countrob, i, GroupName);
  cd(paperfigdir);
  export_fig([name, '_r864.png'], '-r864');
  % export_fig([name, '.pdf']);
end

return
%% Zeichne die Marker in jeweils ein eigenes Bild zum Einfügen in Tabelle
for i = 1:length(RobotGroups) %#ok<UNRCH>
  change_current_figure(1000);clf;
  plot(1,1,RobotGroups{i,4});
  figure_format_publication(gca);
  set(gcf, 'color', 'none'); % transparent. Funktioniert nicht in pdf.
  set(gca, 'Box', 'off');
  set(gca, 'XTICK', [], 'YTICK', []);
  set(get(gca, 'XAXIS'), 'visible', 'off');
  set(get(gca, 'YAXIS'), 'visible', 'off');
  set_size_plot_subplot(1000, ...
    1,1,gca,0,0,0,0,0,0)
  export_fig(fullfile(paperfigdir, sprintf('group%d_marker.pdf', i)));
  export_fig(fullfile(paperfigdir, sprintf('group%d_marker.pdf', i)));
  cd(paperfigdir);
  export_fig(sprintf('group%d_marker.png', i), '-r864');
end
