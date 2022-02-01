% Evaluate Pareto fronts for groups of robots. The groups are selected by
% similar kinematic properties to reduce the amount of solutions.
%
% Preliminaries:
% * run dimensional synthesis with dimsynth_3T1R_example.m
% * Generate robot names with robot_names.m
% * Aggregate results with  eval_figures_pareto.m
% 
% This script is based on the same file eval_figures_pareto_groups.m from 
% https://github.com/SchapplM/robsynth-paper_mhi2021
% (MHI paper "Combined Structural and Dimensional Synthesis of a Parallel
% Robot for Cryogenic Handling Tasks", Schappler et al. 2022)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-01
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear
close all
%% Initialisierung
this_dir = fileparts(which('eval_figures_pareto_groups.m'));
addpath(fullfile(this_dir, '..'));
% Mögliche Zielkriterien: condition, actforce, jointrange, energy
ps = 1; % muss konsistent zu eval_figure_pareto sein
if ps ==  1, pareto_settings = {'chainlength', 'condition'}; end

resdirtotal = ark2022_3T1R_dimsynth_data_dir();
outputdir =this_dir;
datadir = fullfile(fileparts(which('ark2022_3T1R_dimsynth_data_dir.m')),'data');
paperfigdir = fullfile(outputdir, '..', 'paper', 'figures');
% Namen der Roboter laden (enthält auch Eigenschaften zu Kardan-Gelenken)
namestablepath = fullfile(datadir, 'robot_names_latex.csv');
ResTab_NameTrans = readtable(namestablepath, 'Delimiter', ';');
%% Daten laden (als mat)
tmp = load(fullfile(datadir, 'results_all_reps_pareto.mat'));
ResTab_ges = tmp.ResTab_ges;
% ResTab_ges = readtable(fullfile(datadir, 'results_all_reps_pareto.csv'), ...
%   'Delimiter', ';', 'ReadVariableNames', true);

%% Roboter gruppieren
Robots = unique(ResTab_ges.Name);
% Erzeuge die Namen der Gruppen durch weglassen der P- und G-Nummer
RobotGroups = cell(length(Robots),4);
for i = 1:length(Robots)
  [~, ~, ~, ~, ~, ~, ~, PName_Kin, PName_Legs] = parroblib_load_robot(Robots{i});
  RobotGroups{i,1} = [PName_Legs, 'G'];
  % Finde alle Roboter, die zur Gruppe passen und schreibe sie in die
  % zweite Spalte. Marker 'G' am Ende, damit nicht Varianten vereint werden
  RobotGroups{i,2} = Robots(contains(Robots,RobotGroups{i,1}));
  % Dritte Spalte, Marker, ob i.O.-Ergebnisse in Maßsynthese
  RobotGroups{i,3} = 0;
  % Vierte Spalte, Plot-Marker (in allen Bildern konsistent)
  RobotGroups{i,4} = '';
end
[~,I] = unique(RobotGroups(:,1));
RobotGroups = RobotGroups(I,:);
%% Alle Roboter einzeln durchgehen
markerlist = {'x', 's', 'v', '^', '*', '+', '<', '>', 'o', 'd', 'p', 'h'};
colorlist =  {'r', 'g', 'b', 'c', 'm', 'k'};

I_robleg = false(length(Robots), 1);
% Betrachte nur Optimierungsläufe mit i.O.-Ergebnis
I_iO = ResTab_ges.Fval_Opt < 1e3;
% Betrachte nur Optimierungsläufe mit passender Zielfunktion
I_objmatch = contains(ResTab_ges.Zielfunktion,pareto_settings{1}) & ...
             contains(ResTab_ges.Zielfunktion,pareto_settings{2});
change_current_figure(10);clf; hold all;
leghdl = [];
legstr = {};
countrob = 0;
ic = 0;
im = 1;
for i = 1:size(RobotGroups,1)
  GroupName = RobotGroups{i,1};
  I_groupmatch = contains(ResTab_ges.Name, GroupName);
  % Generiere die Marker-Symbole bereits vor der Prüfung, ob die richtigen
  % Zielfunktionen gewählt sind. Dadurch wird sichergestellt, dass die
  % Marker in allen Pareto-Diagrammen gleich sind.
  if any(I_iO&I_groupmatch)
%     im = im+1; % Index für Farben und Marker generieren
    ic = ic+1;
    marker = markerlist{im};
    color = colorlist{ic};
    RobotGroups{i,4} = [marker,color];
    if ic == length(colorlist)
      im = im + 1;
      ic = 0;
    end
  end
  II_Robi = find(I_iO&I_objmatch&I_groupmatch);
  if isempty(II_Robi)
    continue
  end
  fprintf('Gruppe %d/%d (%s): Lade Daten (%d i.O.-Optimierungs-Läufe)\n', i, ...
    size(RobotGroups,1), GroupName, length(II_Robi));
  % disp(ResTab_ges.OptName(II_Robi));
  numrep_i = 0;
  pt_i = cell2table(cell(0,6), 'VariableNames', {'OptName', 'RobName', 'LfdNr', ...
    'ParetoIndNr', 'ChainLength', 'Condition'});
  %% Stelle Pareto-Front aus verschiedenen Durchläufen zusammen
  pf_data = []; % Pareto-Front mit physikalischen Daten. Spalten bezogen auf pareto_settings
  for j = 1:length(II_Robi) % Verschiedene Gut-Durchläufe durchgehen
    OptName = ResTab_ges.OptName{II_Robi(j)};
    RobName = ResTab_ges.Name{II_Robi(j)};
    LfdNr = ResTab_ges.LfdNr(II_Robi(j));
    % Lade die Ergebnisse
    setfile = dir(fullfile(resdirtotal, OptName, '*settings.mat'));
    d1 = load(fullfile(resdirtotal, OptName, setfile(1).name));
    Set_j = d1.Set;
    resfile = fullfile(resdirtotal, OptName, sprintf('Rob%d_%s_Endergebnis.mat', LfdNr, RobName));
    tmp = load(resfile);
    RobotOptRes_j = tmp.RobotOptRes;
    
    % Stelle alle Daten zur Pareto-Front zusammen. Vermeide NaN-Werte, die
    % auftreten, wenn die Nebenbedingungen verletzt werden und dies trotz-
    % dem Teil der Pareto-Front wird.
    I_notnan = all(~isnan(RobotOptRes_j.physval_pareto),2);
    fval_pareto = RobotOptRes_j.fval_pareto(I_notnan,:);
    p_val_pareto = RobotOptRes_j.p_val_pareto(I_notnan,:);
    I_pareto = find(I_notnan); % Indizes der gewählten Partikel in den ursprünglichen Pareto-Fronten
    physval_pareto = RobotOptRes_j.physval_pareto(I_notnan,:);
    % Lade die Detail-Ergebnisse mit allen Zwischenwerten
    resfile2 = fullfile(resdirtotal, OptName, sprintf('Rob%d_%s_Details.mat', LfdNr, RobName));
    if exist(resfile2, 'file')
      t1=tic();
      tmp2 = load(resfile2);
      PSO_Detail_Data_j = tmp2.PSO_Detail_Data;
    else
      % Dummy-Variable anlegen (für Schleife im nächsten Schritt)
      PSO_Detail_Data_j = struct('pval', NaN(3,size(p_val_pareto,2),2), ...
        'physval', NaN(3,size(physval_pareto,2),2), ...
        'fval', NaN(3,size(fval_pareto,2),2));
    end
    
    kk1 = find(strcmp(Set_j.optimization.objective, pareto_settings{1}));
    kk2 = find(strcmp(Set_j.optimization.objective, pareto_settings{2}));
    % Wähle nur Durchläufe, bei denen nur die gewählten Kriterien für
    % die Pareto-Front benutzt wurden. Sonst eher Streudiagramm.
    if isempty(kk1) || isempty(kk2)
      continue % nicht die richtigen Zielkriterien
    end
    
    % Index der Roboter, die Isomorphismen des aktuellen darstellen.
    % Nur sinnvoll, wenn aktueller Roboter ein allgemeiner Typ ist.
    Ir_this = strcmp(ResTab_NameTrans.PKM_Name, RobName); % Index in Namens-Tabelle
    if ~any(Ir_this)
      warning(['Kein Robotername für %s gefunden. Aufruf von robot_names.m ', ...
        'unvollständig.'], RobName);
      continue
    end
    ChainName = ResTab_NameTrans.Chain_Name{Ir_this};
    % Lade Daten der Beinkette
    serroblibpath=fileparts(which('serroblib_path_init.m'));
    NLJ = str2double(ChainName(2));
    mdllistfile_Ndof = fullfile(serroblibpath, sprintf('mdl_%ddof', ...
      NLJ), sprintf('S%d_list.mat',NLJ));
    l = load(mdllistfile_Ndof, 'Names_Ndof', 'AdditionalInfo');
    Ir_db = find(strcmp(l.Names_Ndof, ChainName)); % Index dieses Roboters in Datenbank
    Id_db = find(l.AdditionalInfo(:,3)==Ir_db & l.AdditionalInfo(:,2)); % Index daraus abgeleiteter Varianten
    if ~isempty(Id_db) % es gibt aus diesem Roboter abgeleitete Varianten (und es ist nicht dieser Roboter selbst)
      Chain_DepVarMdl = l.Names_Ndof{Id_db}; % Name der Beinkette der Varianten
      % Index der Roboter, die Isomorphismen des aktuellen darstellen
      Ir_jointiso = find(strcmp(ResTab_NameTrans.Chain_Name, Chain_DepVarMdl) & ...
        ResTab_NameTrans.Gnum==ResTab_NameTrans.Gnum(Ir_this) & ...
        ResTab_NameTrans.Pnum==ResTab_NameTrans.Pnum(Ir_this));
      if isempty(Ir_jointiso), Ir_jointiso = []; end
    else
      Ir_jointiso = [];
    end
    % Merke Indizes von Gelenken, deren DH-Parameter für diesen allgemeinen
    % Roboter nicht Null sein müssen. Wenn sie Null sind, entspricht dieser
    % Roboter der Variante und beide wären Isomorphismen.
    Ijoints_notnull = false(1,NLJ);
    % Alle Gelenke dieses Roboters (bzw. der Beinkette) durchgehen und
    % prüfen, ob daraus abgeleitete Varianten in der Optimierung
    % existieren)
    for II_jointiso = Ir_jointiso % Es kann mehrere Varianten geben
      % Bestimme Zeichenkette der technischen Gelenke (z.B. "PRUR")
      TJ_this = ResTab_NameTrans.Chain_ShortName{Ir_this};
      TJ_jointiso = ResTab_NameTrans.Chain_ShortName{II_jointiso};
      % Roboter-Klasse definieren und deren Methoden nutzen
      RS_this = serroblib_create_robot_class(ChainName);
      RS_var = serroblib_create_robot_class(Chain_DepVarMdl);
      RS_this.set_techjoints(TJ_this);
      RS_var.set_techjoints(TJ_jointiso);
      % Unterschied der technischen Gelenke
      Ijointdiff = RS_this.DesPar.joint_type ~= RS_var.DesPar.joint_type;
      % Bei einer Kardan-Gruppe wird immer der zweite a/d-Parameter zu Null
      % gesetzt (der erste Parameter führt zum Ort des Gelenks)
      ignorenext = false; % Merker zum Finden des jeweils zweiten Eintrags
      for iii = 1:length(RS_var.DesPar.joint_type)
        if ~Ijointdiff(iii), continue; end
        if RS_var.DesPar.joint_type(iii) == 2 && ~ignorenext
          ignorenext = true;
        elseif RS_var.DesPar.joint_type(iii) == 2 && ignorenext
          Ijoints_notnull(iii) = true;
          ignorenext = false;
        end
      end
    end
    % Markiere jeweils die zweiten DH-Parameter einer U-Gelenk-Gruppe
    Ip_DH_Ujoint = false(length(RobotOptRes_j.Structure.varnames),1);
    for ii = 1:length(Ijoints_notnull)
      if Ijoints_notnull(ii)
        Ip_DH_Ujoint(contains(RobotOptRes_j.Structure.varnames, sprintf('d%d', ii))) = true;
        Ip_DH_Ujoint(contains(RobotOptRes_j.Structure.varnames, sprintf('a%d', ii))) = true;
      end
    end
    if any(Ip_DH_Ujoint)
      fprintf(['Gelenkkette %s (%s): Es existiert die Variante %s (%s). Folgende DH-', ...
        'Parameter müssen ungleich Null sein, damit von Variante verschieden: %s\n'], ...
        ChainName, TJ_this, Chain_DepVarMdl, TJ_jointiso, ...
        disp_array(RobotOptRes_j.Structure.varnames(Ip_DH_Ujoint), '[%s]'));
    end
    % Index für Skalierungsparameter
    Ip_scale = strcmp(RobotOptRes_j.Structure.varnames, 'scale');
    % Füge alle temporären Ergebnisse ebenso hinzu. Diese sind nicht
    % Paretodominierend und werden sowieso entfernt.
    % Dadurch ist es möglich, nochmal zu filtern und die Pareto-Front neu
    % aufzustellen.
    for k = 1:size(PSO_Detail_Data_j.pval, 3)
      I_notnan = all(~isnan(PSO_Detail_Data_j.physval(:,:,k)),2);
      % Neue Menge an Partikeln erzeugen (diese Generation dazu).
      fval_pareto = [fval_pareto; PSO_Detail_Data_j.fval(I_notnan,:,k)]; %#ok<AGROW>
      p_val_pareto = [p_val_pareto; PSO_Detail_Data_j.pval(I_notnan,:,k)]; %#ok<AGROW>
      physval_pareto = [physval_pareto; PSO_Detail_Data_j.physval(I_notnan,:,k)]; %#ok<AGROW>
      I_pareto = [I_pareto; find(I_notnan)]; %#ok<AGROW>
      I_select = true(size(p_val_pareto,1),1);
      % Bei Strukturen mit Drehgelenken: Prüfe deren Abstand. Wenn es die
      % gleiche Struktur auch mit Kardan-Gelenken gibt, sind die beiden
      % gleichwertig und der Abstand muss bei dieser allgemeinen Struktur
      % mindestens 20mm betragen (willkürlicher Wert)
      if any(Ip_DH_Ujoint)
        DH_ad_pareto_k_rel = p_val_pareto(:, Ip_DH_Ujoint);
        DH_ad_pareto_k = DH_ad_pareto_k_rel .* ...
          repmat(p_val_pareto(:,Ip_scale), 1, sum(Ip_DH_Ujoint));
        I_select = I_select & all(abs(DH_ad_pareto_k) > 20e-3,2);
      end
      
      % Wende Filter an. Die ersten Generationen werden mehrfach gefiltert.
      % Ist rechentechnisch egal und Code so kompakter.
      p_val_pareto = p_val_pareto(I_select,:);
      physval_pareto = physval_pareto(I_select,:);
      I_pareto = I_pareto(I_select);
      fval_pareto = fval_pareto(I_select,:);
      % Bereits hier wieder Pareto-Dominanz prüfen. Sonst führt das weiter
      % unten zur Überlastung durch die Vielzahl an Partikeln
      if size(physval_pareto,1) > 1
        Idom_ges = pareto_dominance(physval_pareto);
        physval_pareto = physval_pareto(~Idom_ges,:);
        p_val_pareto = p_val_pareto(~Idom_ges,:);
        fval_pareto = fval_pareto(~Idom_ges,:);
        I_pareto = I_pareto(~Idom_ges);
      end
    end
    assert(size(physval_pareto,1)==size(p_val_pareto,1), 'Falsche Dimension');
    assert(size(physval_pareto,1)==size(fval_pareto,1), 'Falsche Dimension');
    assert(size(physval_pareto,1)==size(fval_pareto,1), 'Falsche Dimension');
    pf_data = [pf_data; physval_pareto(:,[kk1,kk2])]; %#ok<AGROW>
    row_i = cell(size(physval_pareto,1),6);
    row_i(:,1:3) = repmat({OptName,RobName,LfdNr},size(row_i,1),1);
    for k = 1:size(p_val_pareto,1)
      row_i{k,4} = I_pareto(k); % Pareto Index. TODO: Bezug stimmt aktuell nicht
      row_i{k,5} = physval_pareto(k,kk1); % Kriterium 1: Chainlength
      row_i{k,6} = physval_pareto(k,kk2); % Kriterium 2: Condition
      assert(all(~isnan(physval_pareto(k,:))), 'physval_pareto ist NaN.');
    end
    pt_i = [pt_i; row_i]; %#ok<AGROW>
    numrep_i = numrep_i + 1;
  end
  if ~isempty(pf_data)
    [~, Ikk] = sort(pf_data(:,1)); % Sortiere nach erstem Kriterium
    pt_i = pt_i(Ikk,:);
    pf_data = pf_data(Ikk,:);
  end
  % Erstelle eine einzige Pareto-Front. Die Fronten mehrere Durchläufe sind
  % durch die heuristische Optimierung gegenseitig dominant.
  % Definition Pareto-Front: Siehe CoelloPulLec2004 Gl. (1)-(6)
  Idom_ges = pareto_dominance(pf_data);
  fprintf(['Durch mehrfache Durchführung der Optimierung müssen %d/%d Partikel ', ...
    'von der Pareto-Front entfernt werden.\n'], sum(Idom_ges), length(Idom_ges));
  pf_data = pf_data(~Idom_ges,:);

  % Markiere diese Gruppe als i.O. (da Ergebnisse vorliegen)
  RobotGroups{i,3} = size(pf_data,1);

  pt_i = pt_i(~Idom_ges,:);
  % Speichere die Ergebnisse der Daten für diesen Roboter (bzw. die Gruppe)
  writetable(pt_i, fullfile(datadir, ...
    sprintf('group_%s_paretofront.csv', GroupName)), 'Delimiter', ';');
  save(fullfile(datadir, sprintf('group_%s_paretofront.mat', GroupName)), 'pt_i');
  
  if isempty(pf_data)
    % Durch Filterkriterien wird die Gruppe doch wieder aussortiert.
    continue;
  end
  % Ab hier ist ein Roboter erfolgreich
  countrob = countrob + 1; % Für Legende: Nur erfolgreiche PKM zählen
  I_robleg(i) = true;
  %% Zeichne die Ergebnisse in das Bild ein
  % Bild mit physikalischen Werten
  
  % Zeichne eine Linie, falls mehrere Punkte beieinander sind. Ansonsten
  % zeichne einzelne Marker
  % Abstand zwischen zwei Punkten, bei dem die Linie unterbrochen werden
  % soll. Wird per Hand eingestellt. Annahme: Dann kann man nicht davon
  % ausgehen, dass kleine Änderungen der Kinematikparameter möglich sind.
  % (zusammenhängender Lösungsraum)
  % Zahlen beziehen sich auf die Werte vor der Einheitenkorrektur (deg, µm)
  if ps == 1
    maxgaplength_x = 0.5; % bezogen auf x-Achse, in m. TODO!
  else
    error('Value for ps not set yet');
  end
  xminmax = minmax2(pf_data(:,1)');
  if ps == 1
    dx = 5e-3; % Prüfe in 5mm-Schritten. TODO!
  else
    error('Value for ps not set yet');
  end
  gapdata_x = (xminmax(1):dx:xminmax(2))';
  gapdata_y = NaN(length(gapdata_x),1);
  for jj = 1:length(gapdata_x)
    mindist_jj_x = min(abs(gapdata_x(jj) - pf_data(:,1)));
    if mindist_jj_x > maxgaplength_x/2 % der Abstand gilt immer zu beiden Seiten
      % Setze NaN, damit hier eine Lücke gemacht wird
      gapdata_y(jj) = NaN;
    else
      gapdata_y(jj) = 0; % Setze Null, damit Wert am Ende entfernt wird
    end
  end
  pf_plotdata = [pf_data; [gapdata_x(gapdata_y~=0), gapdata_y(gapdata_y~=0)]];
  % Daten neu stetig sortieren
  [~,II] = sort(pf_plotdata(:,1));
  pf_plotdata = pf_plotdata(II,:);
  % Dünne die Daten aus, damit die Marker sich nicht komplett gegenseitig
  % überdecken (nur auf der Linie)
  % minimaler Abstand zwischen zwei Markern (x-, y-Achse). Muss per Hand so
  % eingestellt werden, dass es mit der Skalierung des Bildes passt.
  % Zahlen beziehen sich auf die Werte vor der Einheitenkorrektur (deg, µm)
  if ps == 1
    minmarkerdist = [0.3,5]; % je größer der Wert desto enger die Marker. TODO
  else
    error('Value for ps not set yet');
  end
  I_pm = true(size(pf_plotdata,1),1); % true, falls Marker für Partikel gezeichnet wird
  for jj = 1:length(I_pm)
    if jj == 1 || isnan(pf_plotdata(jj-1,2)) || jj~=length(I_pm) && isnan(pf_plotdata(jj+1,2))
      continue % lasse die 1, zeichne Marker
    end
    I_left = (1:length(I_pm))'<jj;
    I_nan = isnan(pf_plotdata(:,2)); % NaN-Marker
    I_leftnan = [I_nan(2:end);false]; % Eine Position links von Lücke (entspricht rechtem Marker einer Linie
    I_rightnan = [false; I_nan(1:end-1)];
    % Prüfe den Abstand des aktuellen Markers zu allen anderen
    markerdist_jj_x = min(abs(pf_plotdata(jj,1)-pf_plotdata(I_pm&(I_left),1)));
    markerdist_jj_y = min(abs(pf_plotdata(jj,2)-pf_plotdata(I_pm&(I_left|I_leftnan|I_rightnan),2)));
    if markerdist_jj_x < minmarkerdist(1) || markerdist_jj_y < minmarkerdist(2)
      % Dieser Marker ist zu nah an einem anderen dran. Zeichne ihn nicht.
      I_pm(jj) = false;
    end
  end
  % Linienstil bestimmen und zeichnen
  if all(isnan(gapdata_y)) % Alle Werte sind einzeln, keine Linie zu zeichnen
    linestyle = '';
  else % TODO: Hier kann noch eine andere Strichelung genommen werden
    linestyle = '-';
  end
  [obj_units, objscale] = cds_objective_plotdetails(Set_j);
  obj_units{strcmp(Set_j.optimization.objective, 'chainlength')} = 'm';
  objscale(strcmp(Set_j.optimization.objective, 'chainlength')) = 1;
  % Zeichne durchgezogene Linie
  plot(objscale(1)*pf_plotdata(:,1), objscale(2)*pf_plotdata(:,2), [color,linestyle]);
  % Zeichne Marker
  plot(objscale(1)*pf_plotdata(I_pm,1), objscale(2)*pf_plotdata(I_pm,2), [color,marker]);
  % Dummy-Handle für Legende
  hdl = plot(0, NaN, [color, linestyle,marker]);
  xlabel(sprintf('%s in %s', pareto_settings{1}, obj_units{1}));
  ylabel(sprintf('%s in %s', pareto_settings{2}, obj_units{2}));
  grid on;
  leghdl(countrob) = hdl; %#ok<SAGROW>
  legstr{countrob} = sprintf('%s; %d Wdh.', GroupName, numrep_i); %#ok<SAGROW>
end
set(10, 'numbertitle', 'off');
title(sprintf('Pareto-Front %s vs %s (Fitness-Werte physikalisch)', pareto_settings{1}, pareto_settings{2}));
lh = legend(leghdl, legstr);
set(10, 'name', sprintf('pareto_groups_%s_%s', pareto_settings{1}(1:4), pareto_settings{2}(1:4)));
saveas(10, fullfile(outputdir, sprintf('figure_pareto_groups_%s_%s.fig', pareto_settings{1}, pareto_settings{2})));
export_fig(10, fullfile(outputdir, sprintf('figure_pareto_groups_%s_%s.pdf', pareto_settings{1}, pareto_settings{2})));

%% Gruppen abspeichern
save(fullfile(datadir, 'robot_groups.mat'), 'RobotGroups');
maxgroupsize = max(cellfun(@length,RobotGroups(:,2)));
varnames = cell(1,maxgroupsize);
varnames{1} = 'Number';
varnames{2} = 'GroupName';
varnames{3} = 'Chain_ShortName';
for i = 1:maxgroupsize
  varnames{3+i} = sprintf('RobName%d', i);
end
GroupTab = cell2table(cell(0,3+maxgroupsize), 'VariableNames', varnames);
for i = 1:size(RobotGroups,1)
  grrow = cell(1,3+maxgroupsize);
  grrow{1} = i;
  grrow{2} = RobotGroups{i,1};
  % Suche den verkürzten Namen der kinematischen Kette
  I = strcmp(ResTab_NameTrans.PKM_Name, RobotGroups{i,2}{1});
  if any(I)
    grrow(3) = ResTab_NameTrans.Chain_ShortName(I);
    for j = 1:length(RobotGroups{i,2})
      grrow{3+j} = RobotGroups{i,2}{1};
    end
  else
    warning('Group %d (%s) is empty', i, RobotGroups{i,1});
  end
  GroupTab = [GroupTab; grrow]; %#ok<AGROW>
end
writetable(GroupTab, fullfile(datadir, 'robot_groups.csv'), 'Delimiter', ';');

%% Für Paper formatieren
% uiopen(fullfile(outputdir, sprintf('figure_pareto_groups_%s_%s.fig', ...
%   pareto_settings{1}, pareto_settings{2})), 1);

% Legende und Roboternamen aktualisieren
legstr2 = cell(length(legstr),1);
% Ersetze den ursprünglichen Legendennamen durch die Bezeichnung im Paper
for i = 1:length(legstr)
  for j = 1:size(RobotGroups,1) % muss nicht identisch mit Legende sein
    if contains(legstr{i}, RobotGroups{j,1})
      % Suche Latex-Namen
      for k = 1:length(RobotGroups{j,2})
        I_k = strcmp(ResTab_NameTrans.PKM_Name, RobotGroups{j,2}{k});
        legstr2{i} = sprintf('%d: 4-%s (%s)', i, ResTab_NameTrans.Chain_Structure_Act{I_k}, ...
          ResTab_NameTrans.Chain_ShortName{I_k});
        break;
      end
      break;
    end
    if ~isempty(legstr2{i})
      break;
    end
  end
end

if ps == 1
  xlabel('Sum of leg chain lengths in m');
else
  error('Value for ps not set yet');
end
if ps == 1
  ylabel('Jacobian condition number');
else
  error('Value for ps not set yet');
end
title('');

if ps == 1 % Manuell einstellen
%   axis auto
  set(gca, 'yscale', 'log');
  xlim([8, 22]);
  ylim([3, 1e4]);
else
  error('Value for ps not set yet');
end
grid on
figure_format_publication(gca);
set_size_plot_subplot(10, ...
  11.6,4.5,gca,...
  0.08,0.005,0.04,0.10,0,0)

lh = legend(leghdl, legstr2, 'interpreter', 'latex');
% Move legend manually and obtain position by `get(lh, 'position')`. Then
% set these values below
set(lh, 'Position', [0.6330    0.1417    0.3472    0.7882]);

% Ziehe die x-Beschriftung etwas nach oben und entferne dort die Ticklabel
xtl = get(gca, 'xticklabel');
xtl(3:4) = {''};
set(gca, 'xticklabel', xtl);
xlhdl = get(gca, 'xlabel');
[X_off, X_slope] = get_relative_position_in_axes(gca, 'x');
[Y_off, Y_slope] = get_relative_position_in_axes(gca, 'y');
xlpos = get(xlhdl, 'Position');
if ps == 1
  xlpos(1) = X_off+X_slope*(-0.4);
else
  error('Value for ps not set yet');
end
xlpos(2) = Y_off*Y_slope^(-1/2+0.5-0.02);
set(xlhdl, 'Position', xlpos);

% Speichern und Abschluss
name = sprintf('paperfigure_pareto_groups_%s_%s', pareto_settings{1}, pareto_settings{2});
saveas(10, fullfile(outputdir, sprintf('%s.fig', name)));
export_fig(10, fullfile(outputdir, sprintf('%s.pdf', name)));
export_fig(10, fullfile(paperfigdir, 'pareto_all.pdf'));
fprintf('Auswertung %s gespeichert\n', name);
% Vergrößerung zusätzlich speichern. In einem Bereich ist es zu
% teilweise zu unübersichtlich (nicht mehr genutzt).
return
delete(lh);
xlabel(''); ylabel('');
if ps == 11
  xlim([26 29]);
  ylim([38 50]);
else
  error('Value for ps not set yet');
end
set_size_plot_subplot(10, ...
  4,4,gca,...
  0.10,0.02,0.02,0.12,0,0)
name = sprintf('paperfigure_pareto_groups_%s_%s_detail', pareto_settings{1}, pareto_settings{2});
saveas(10, fullfile(outputdir, sprintf('%s.fig', name)));
export_fig(10, fullfile(outputdir, sprintf('%s.pdf', name)));
export_fig(10, fullfile(paperfigdir, 'pareto_all_detail.pdf'));


