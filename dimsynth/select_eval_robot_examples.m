% Select robots for a detailed inspection in the paper
% (mainly detail plots)
% influences robot_images.m 
% take robot from the Pareto front and re-compute the fitness function
% 
% This script is based on the same file select_eval_robot_examples.m from 
% https://github.com/SchapplM/robsynth-paper_mhi2021
% (MHI paper "Combined Structural and Dimensional Synthesis of a Parallel
% Robot for Cryogenic Handling Tasks", Schappler et al. 2022)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-01
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear
%% Benutzereingaben
recalc_fitnessfcn = false; % Neuberechnung der Fitness-Funktion (optional)
% % Wähle das Partikel, dessen Positionsfehler am nächsten am gewünschten
% % Wert ist. Wähle nicht Roboter mit möglichst kleiner Antriebskraft, da
% % dann die Beinketten sehr lang werden und die Roboter komisch aussehen.
% posacc_sel = 40e-6;
% Compute redundancy map and save it next to the results.
recreate_redundancy_map = false;
regenerate_templates = false; %#ok<*UNRCH> % nur bei erstem Aufruf notwendig.
%% Sonstige Initialisierung
if isempty(which('ark2022_3T1R_dimsynth_data_dir'))
  error(['You have to create a file ark2022_3T1R_dimsynth_data_dir pointing to the ', ...
    'directory containing the results of the dimensional synthesis']);
end
importdir = ark2022_3T1R_dimsynth_data_dir();
datadir = fullfile(fileparts(which('ark2022_3T1R_dimsynth_data_dir.m')),'data');
tmp = load(fullfile(datadir, 'results_all_reps_pareto.mat'));
ResTab = tmp.ResTab_ges;
tmp = load(fullfile(datadir, 'robot_groups.mat'));
RobotGroups = tmp.RobotGroups;

% Debug: Gruppen manuell finden
for kk = 1:size(RobotGroups)
  if any(strcmp(RobotGroups{kk,2}, 'P4PRRRR6G2P2A1'))
    fprintf('Gruppe %d: %s\n', kk, RobotGroups{kk,1})
    break;
  end
end

%% Alle Gruppen durchgehen
for i = 1:size(RobotGroups,1)
  GroupName = RobotGroups{i,1};
  if RobotGroups{i,3} == 0, continue; end % keine Ergebnisse vorliegend
  fprintf('Lade Daten für PKM-Gruppe %d/%d (%s)\n', i, size(RobotGroups,1), GroupName);
  data_i = load(fullfile(datadir, sprintf('group_%s_paretofront.mat', GroupName)));
  % Manuelle Anpassung der zu findenden Roboter
  II = (1:length(data_i.pt_i.RobName))';
  % Suche den Eintrag mit der kürzesten Beinkettenlänge bei gleichzeitig
  % niedriger Konditionszahl (Grenzwert 200 manuell gewählt)
  I_lowcond = data_i.pt_i.Condition < 200;
  if any(I_lowcond)
    II_lc = II(I_lowcond);
    [~,iinearest] = min(abs(data_i.pt_i.ChainLength(II_lc)));
    inearest = II_lc(iinearest);
  else
    % Kein Roboter hat eine akzeptable Konditionszahl. Nehme direkt den
    % Roboter mit der kürzesten kinematischen Kette
    [~,iinearest] = min(abs(data_i.pt_i.ChainLength(II)));
    inearest = II(iinearest);
  end
  % Lade Daten für diesen Roboter aus den Ergebnissen
  Ipar = data_i.pt_i.ParetoIndNr(inearest);
  OptName = data_i.pt_i.OptName{inearest};
  RobName = data_i.pt_i.RobName{inearest};
  LfdNr = data_i.pt_i.LfdNr(inearest);
  fprintf('Wähle Opt. %s, Rob. %d, %s, Partikel %d\n', OptName, LfdNr, RobName, Ipar);
  setfile = dir(fullfile(importdir, OptName, '*settings.mat'));
  d1 = load(fullfile(importdir, OptName, setfile(1).name));
  Set_i = cds_settings_update(d1.Set);
  resfile = fullfile(importdir, OptName, sprintf('Rob%d_%s_Endergebnis.mat', LfdNr, RobName));
  tmp = load(resfile);
  RobotOptRes_i = tmp.RobotOptRes;
  resfile2 = fullfile(importdir, OptName, sprintf('Rob%d_%s_Details.mat', LfdNr, RobName));
  if exist(resfile2,'file')
    % Lade auch die Detail-Ergebnisse. Geht nicht ohne, wenn die
    % Detail-Ergebnisse genutzt wurden, um die Pareto-Front anzupassen.
    % (z.B. durch nachträgliche Filterung nach zusätzlichen Kriterien
    tmp = load(resfile2);
    PSO_Detail_Data_i = tmp.PSO_Detail_Data;
    RobotOptDetails_i = tmp.RobotOptDetails;
  else
    PSO_Detail_Data_i = [];
    RobotOptDetails_i = [];
  end
  % Nummer der Zielfunktion "Positionsfehler"
  kk1 = find(strcmp(Set_i.optimization.objective, 'chainlength'));
  kk2 = find(strcmp(Set_i.optimization.objective, 'condition'));
  % Prüfe, ob Werte zueinander passen. Kann nicht mehr direkt die
  % Pareto-Fronten nehmen, da die Pareto-Fronten neu gebildet werden. Es
  % muss immer das passende Partikel in den Zwischenergebnissen gesucht
  % werden.
  if abs(RobotOptRes_i.physval_pareto(Ipar,kk1) - data_i.pt_i.ChainLength(inearest)) < 1e-10
    % Der Gesuchte Wert liegt ganz normal auf der Pareto-Front aus dem
    % Optimierungsergebnis
    % Lade weitere Daten aus der Ergebnis-Datei (aus Endergebnis-Variable)
    pval = RobotOptRes_i.p_val_pareto(Ipar,:)';
    physval = RobotOptRes_i.physval_pareto(Ipar,:)';
    fval = RobotOptRes_i.fval_pareto(Ipar,:)';
    pval_desopt = RobotOptRes_i.desopt_pval_pareto(Ipar,:)';
    q0 = RobotOptRes_i.q0_pareto(Ipar,:)';
  elseif ~isempty(PSO_Detail_Data_i)
    % Suche in den Detail-Daten. Nehme beide Kriterien, damit die Daten
    % eindeutig sind. Ansonsten können mehrere Parametersätze teilweise die
    % gleichen Zielkriterien erzeugen
    ChainlengthMatrix = reshape(PSO_Detail_Data_i.physval(:,kk1,:), ...
      size(PSO_Detail_Data_i.physval,1), size(PSO_Detail_Data_i.physval,3));
    ConditionMatrix = reshape(PSO_Detail_Data_i.physval(:,kk2,:), ...
      size(PSO_Detail_Data_i.physval,1), size(PSO_Detail_Data_i.physval,3));
    k = find( ...
      abs(ChainlengthMatrix(:)  -data_i.pt_i.ChainLength(inearest))<1e-10 & ...
      abs(ConditionMatrix(:)-data_i.pt_i.Condition(inearest))<1e-10);
    if length(k) > 1 % Es gibt mehrere Partikel mit genau diesen Werten für die Zielkriterien
      % Prüfe, ob wenigstens die Parameter unterschiedlich sind.
      pval_all = NaN(length(k), size(PSO_Detail_Data_i.pval,2));
      for ii = 1:length(k)
        [ii_ind,ii_gen] = ind2sub(fliplr(size(PSO_Detail_Data_i.comptime)),k(ii));
        pval_all(ii,:) = PSO_Detail_Data_i.pval(ii_ind,:,ii_gen);
      end
      if size(unique(pval_all,'rows'), 1) == 1
        k = k(1); % sind unterschiedlich. Also identische Roboter.
      else % nicht unterschiedlich. Prinzipielles Problem (redundante Parameter)
        error('Suchkriterium in Daten ist nicht eindeutig');
      end
    end
    [k_ind,k_gen] = ind2sub(fliplr(size(PSO_Detail_Data_i.comptime)),k);
    physval = PSO_Detail_Data_i.physval(k_ind,:,k_gen)';
    if abs(physval(kk1)-data_i.pt_i.ChainLength(inearest))>1e-10
      error('Gesuchter Wert konnte nicht gefunden werden. Logik-Fehler');
    end
    fval = PSO_Detail_Data_i.fval(k_ind,:,k_gen)';
    pval = PSO_Detail_Data_i.pval(k_ind,:,k_gen)';
    pval_desopt = PSO_Detail_Data_i.desopt_pval(k_ind,:,k_gen)';
    q0 = PSO_Detail_Data_i.q0_ik(k_ind,:,k_gen)';
  else
    error(['Ergebnis-Partikel liegt nicht in der finalen Pareto-Front ', ...
      'der Optimierung. Detail-Auswertung notwendig. Datei aber nicht da: %s'], ...
      resfile2);
  end

  if ~isempty(PSO_Detail_Data_i)
    % Lese die IK-Anfangswerte aus den Ergebnissen aus (sind indirekt in ge-
    % speicherten Zwischenwerten enthalten)
    [k_gen, k_ind] = cds_load_particle_details(PSO_Detail_Data_i, fval);
    pval_test = PSO_Detail_Data_i.pval(k_ind,:,k_gen)' - pval;
    if any(abs(pval_test))
      error('Indizies für Generation/Individuum stimmen nicht (Parameter nicht gleich)');
    end
    
  end
  %% Nachrechnen der Fitness-Funktion
  % um die Gelenkwinkel aus der IK zu erhalten. Annahme: Die gespeicherten 
  % Detail-Informationen stehen aus Speicherplatzgründen nicht zur Verfügung.
  parroblib_addtopath({RobName}); % Für Ausführung der Fitness-Fcn
  if regenerate_templates
    parroblib_create_template_functions({RobName}, false); % Für Erstellung fehlender Dateien
    R_test = parroblib_create_robot_class(RobName, 1, 1);
    R_test.fill_fcn_handles(true, true); % Zur Kompilierung fehlender Funktionen zum Nachrechnen der Fitness-Funktion
  end

  [R, Structure] = cds_dimsynth_robot(Set_i, d1.Traj, d1.Structures{LfdNr}, true);
  % Fitness-Funktion neu definieren (mit weniger Log-Ausgaben)
  Set = Set_i;
  kk1 = strcmp(Set.optimization.objective,'chainlength');
  Set.general.plot_details_in_fitness = 0; % debug: 1e10
  Set.general.save_robot_details_plot_fitness_file_extensions = {};
  Set.general.verbosity = 4; % debug: 3
  Set.optimization.resdir = ark2022_3T1R_dimsynth_data_dir();
  if recreate_redundancy_map
    resdir_tmp = fullfile(Set.optimization.resdir, Set.optimization.optname, ...
      'tmp', sprintf('%d_%s', Structure.Number, Structure.Name));
    mkdirs(resdir_tmp);
    Set.general.debug_taskred_perfmap = 1;
  end
  Set.general.save_robot_details_plot_fitness_file_extensions = {'png', 'fig'};
  % Schreibe temporäre Daten dahin, wo auch die Ursprungsdaten liegen.
  % Überschreibt nichts.

  cds_log(0, '', 'init', Set);
  % IK-Anfangswerte für dieses Partikel setzen
  for iLeg = 1:R.NLEG
    R.Leg(iLeg).qref = q0(R.I1J_LEG(iLeg):R.I2J_LEG(iLeg));
  end
  % Funktionsaufruf siehe cds_check_results_reproducability.m
  Structure_tmp = Structure;
  Structure_tmp.calc_dyn_act = Structure.calc_dyn_act | Structure.calc_dyn_reg;
  Structure_tmp.calc_spring_act = Structure.calc_spring_act | Structure.calc_spring_reg;
  Structure_tmp.calc_spring_reg = false;
  Structure_tmp.calc_dyn_reg = false;
  % Erzwinge Prüfung dieses Anfangswerts für Trajektorie (falls IK anderes
  % Ergebnis hat). Diese Option sollte nicht notwendig sein. Ist sie leider
  % teilweise.
  Structure_tmp.q0_traj = q0;
  clear cds_save_particle_details cds_fitness
  for repro_retry_iter = 1:2 % TODO: Muss eigentlich direkt funktionieren. Warum nicht beim ersten Mal?
    [fval_i_test, physval_i_test, Q] = cds_fitness(R, Set,d1.Traj, ...
      Structure_tmp, pval, pval_desopt);
    PSO_Detail_Data_tmp = cds_save_particle_details(Set, R, 0, 0, NaN, NaN, NaN, NaN, 'output');
    condJ = PSO_Detail_Data_tmp.constraint_obj_val(1, 4, 1);  
    if any(fval_i_test > 1e3)
      warning('Die Nebenbedingungen wurden bei erneuter Prüfung verletzt');
      Set.optimization.pos_ik_tryhard_num = 100; % Viel mehr IK-Neuversuche, damit Ergebnis auf jeden Fall reproduzierbar
      % Versuche nochmal neu, die Fitness-Funktion zu berechnen
%       if repro_retry_iter == 2
%         error(['Keine gültige Gelenkwinkel berechnet. Parallelität der ', ...
%           'Gelenke und damit Name nicht bestimmbar.'])
%       end
    end
  end
  if any(fval_i_test > 1e3)
    continue % Roboter überspringen. TODO: Soll nicht passieren
  end
  q0_neu = PSO_Detail_Data_tmp.q0_ik(1,:,1)';
  test_fval = fval - fval_i_test;
  test_physval = physval - physval_i_test;
  if abs(test_fval(kk1)) > 1e-6 % Durch geänderte Implementierung möglich
    warning(['Die Antriebskraft hat einen anderen Wert beim neu nachrechnen ', ...
      '(%1.1f vs %1.1f; Diff %1.1e). Physikalischer Wert: %1.4f vs %1.4f ', ...
      '(Diff %1.1e). IK hat anderes Ergebnis!'], fval_i_test(kk1), fval(kk1), ...
      test_fval(kk1), physval_i_test(kk1), physval(kk1), test_physval(kk1));
    if abs(test_physval(kk1)) > 5 % 5N wird noch für das Bild toleriert (Annahme: Auf Cluster war es richtig)
      error('Der Fehler ist zu groß. Das lässt sich nicht mehr mit Zufallszahlen-Toleranz erklären');
    end
  end
  if any(abs(test_fval(~kk1)) > 1e-5)
    warning(['Andere Zielfunktionen haben einen anderen Wert beim neu nachrechnen: ', ...
      '[%s] vs [%s]. Physikalisch: [%s] vs [%s]'], ...
      disp_array(fval_i_test(~kk1)','%1.6e'), disp_array(fval(~kk1)','%1.6f'), ...
      disp_array(physval_i_test(~kk1)','%1.3e'), disp_array(physval(~kk1)','%1.3e'));
  end

  %% Abschließende Berechnungen und Abspeichern
  % Speichere die Trajektorie in der Variable X (für späteres Plotten)
  Traj_0 = cds_transform_traj(R, d1.Traj);
  X = Traj_0.X;

  % Ergebnis-Variable abspeichern
  objective_names = Set_i.optimization.objective;
  save(fullfile(datadir, sprintf('detail_result_group_%s.mat', GroupName)), ...
    'R', 'pval', 'fval', 'physval', 'condJ', 'Structure', 'Q', 'X', ... % Daten zum Ergebnis
    'objective_names', 'OptName', 'RobName', 'LfdNr', 'Ipar'); % Herkunft des Ergebnisses
end
fprintf('Je ein Ergebnis aus %d verschiedenen Gruppen ausgewählt und gespeichert\n', size(RobotGroups,1));
