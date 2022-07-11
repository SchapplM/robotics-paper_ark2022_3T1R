% Video zum Zusammenbau einer 3T1R-PKM in der inversen Kinematik.
% Zeige den Einfluss der Redundanzauflösung mit verschiedenen Kriterien

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-06
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

if isempty(which('parroblib_path_init.m'))
  warning('Repo mit parallelen Robotermodellen ist nicht im Pfad. Beispiel nicht ausführbar.');
  return
end

usr_create_animation = true;
usr_debug_plot_ik = true;

%% Roboter-Klasse vorbereiten
% Ergebnisse laden
datadir = fullfile(fileparts(which('ark2022_3T1R_dimsynth_data_dir')), 'data');
respath = fullfile(fileparts(which('ark2022_3T1R_dimsynth_data_dir')), 'presentation');
GroupName = 'P4RRRRR7G'; % 'P4RRRRR5G'; P4RRRRR7G'; P4RRRRR10G';
erg = load(fullfile(datadir, sprintf('detail_result_group_%s.mat', GroupName)));

% Verschiedene IK-Zielkriterien durchgehen
ikobj_all = {'none', 'jac_cond', 'ikjac_cond', 'qlim_par', 'coll_par'};
for i_ikobj = 1:length(ikobj_all)

%% Roboter initialisieren
R = erg.R;
Name = R.mdlname;
% Dicke der Segmente reduzieren
for kk = 1:R.NLEG
  R.Leg(kk).DesPar.seg_par(:,1) = 1e-3;
  R.Leg(kk).DesPar.seg_par(:,2) = 1e-2;
end
parroblib_addtopath({R.mdlname});

X = erg.X(1,:)';

% Start-Gelenkwinkel wählen. TODO: Mit dieser Konfiguration aktuell noch
% falsche Konfiguration der ersten Beinkette als Ergebnis
% q0_leg = zeros(R.Leg(1).NJ,1);
% q0_leg(1) = -pi/4;
% q0_leg(2) = 0;
% q0_leg(2) = 0;
% q0 = repmat(q0_leg, R.NLEG, 1);
% Zufällige Anfangswerte
q0 = erg.Q(1,:)';
rng(0); % Reproduzierbar machen
q0 = q0 + (-0.5+rand(length(q0),1))*2*60*pi/180;
% q0(1:R.Leg(1).NJ:R.NJ) = -pi/2; % nach außen drehen

qlim_pkm = cat(1, R.Leg.qlim);

change_current_figure(2); clf; hold on; grid on; % Bild der Entwurfsparameter
xlabel('x in m');ylabel('y in m');zlabel('z in m'); view(3);
  s_plot = struct( 'ks_legs', [R.I1L_LEG; R.I1L_LEG+1; R.I2L_LEG], ...
    'ks_platform', [1:R.NLEG,R.NLEG+2], 'straight', 1, 'mode', 4);
R.plot( q0, X, s_plot );
title(sprintf('%s in Startkonfiguration', Name));

%% IK zum Startpunkt lösen
t1 = tic();
s = struct('retry_limit', 0, 'normalize', false);
% Keine Gelenkgrenzen in Iterationen beachten. Führt zu stufenweisem
% Verlauf
% s.scale_lim = 0.7;
s.wn = zeros(R.idx_ik_length.wnpos, 1);
% Zielkriterium definieren
if strcmp(ikobj_all{i_ikobj}, 'none')
  % Auslassen
else
  s.wn(R.idx_ikpos_wn.(ikobj_all{i_ikobj})) = 1;
end
s.maxrelstep = 0.01; % sehr feinschrittige Bewegungen (für flüssige Animation)
s.n_max = 5e3;
[q, Phi, ~, Stats] = R.invkin4(X, q0, s);
fprintf('%s: IK berechnet. Kriterium %s. Dauer %1.1fs\n', R.mdlname, ikobj_all{i_ikobj}, toc(t1));

% Berechne die Jacobi-Matrix bezogen auf 3T3R
X_neu = R.fkineEE_traj(q')'; % damit PKM-Jacobi mit Methode 4 stimmt (Orientierung neu bei 3T2R)

% Extrahiere Gelenkverlauf aus IK-Statistik
Q_t_Anfahrt = Stats.Q;
Q_t_norm = (Stats.Q-repmat(qlim_pkm(:,1)',size(Stats.Q,1), 1)) ./ ...
            repmat(qlim_pkm(:,2)'-qlim_pkm(:,1)',size(Stats.Q,1), 1);
Q_t_Anfahrt = Q_t_Anfahrt(1:Stats.iter+1,:);
assert(all(abs(Q_t_Anfahrt(end,:)-q')<1e-3), ...
  'IK-Ergebnis passt nicht zu Anfahrt-Trajektorie');

% Alternative 1: Die Plattform-Pose basierend auf der ersten Beinkette
% rechnen. Führt dazu, dass die Plattform an Beinkette 1 hängt. Dann nicht
% so gut sichtbar, wie die Plattform sich dreht
% X_t_Anfahrt = R.fkineEE_traj(Stats.Q(1:Stats.iter+1,:));
% Alternative 2: Plattform steht an Zielort und dreht sich nur bei der
% Animation.
X_t_Anfahrt = repmat(X_neu, size(Q_t_Anfahrt,1), 1);
fprintf('EE-Drehung nach Ende der IK: %1.1f deg\n', 180/pi*X_t_Anfahrt(end,6));

% Verlauf der IK-Konvergenz prüfen
if usr_debug_plot_ik
  fighdl_diag = change_current_figure(10);clf;
  subplot(2,3,1);
  plot(Q_t_norm(1:max(Stats.iter),:));
  ylabel('Q norm'); grid on;
  subplot(2,3,2);
  plot(Stats.PHI(1:max(Stats.iter),:));
  ylabel('PHI'); grid on;
  subplot(2,3,3);
  plot(X_t_Anfahrt(1:max(Stats.iter),1:3));
  ylabel('x (transl)'); grid on;
  subplot(2,3,4);
  plot(X_t_Anfahrt(1:max(Stats.iter),4:6));
  ylabel('x (rot)'); grid on;
  subplot(2,3,5); hold on
  plot(Stats.condJ(1:max(Stats.iter)));
  if isfield(Stats, 'h')
    plot(Stats.h(1:max(Stats.iter),5));
    legend({'Phi_q', 'J_{PKM}'});
  end
  ylabel('cond(J)'); grid on;
  subplot(2,3,6);
  if isfield(Stats, 'h')
    plot(Stats.h(1:max(Stats.iter),:));
    ylabel('h'); grid on;
  end
  linkxaxes
  sgtitle(sprintf('IK-Diagnose'));
  saveas(fighdl_diag, fullfile(respath, sprintf( ...
    '%s_IK_Stats_Opt_%s.fig',Name, ikobj_all{i_ikobj})));
  change_current_figure(3); clf; hold on; grid on; % Bild der Entwurfsparameter
  xlabel('x in m'); ylabel('y in m'); zlabel('z in m'); view(3);
  s_plot = struct( 'ks_legs', [R.I1L_LEG; R.I2L_LEG], ...
    'straight', 1, 'mode', 4, 'ks_platform', [1:R.NLEG,R.NLEG+2]);
  R.plot( q, X_neu, s_plot );
  title(sprintf('%s nach IK', Name));
end
if any(Stats.iter == s.n_max)
  warning(['Es wurden %d Iterationen gebraucht. Starkes Indiz für ', ...
    'Oszillationen am Ende oder schlechte Konvergenz.'], max(Stats.iter));
end
if any(abs(Phi)>1e-8)
  warning('IK nicht erfolgreich');
  continue % Kein Video hierfür zeichnen (aber kein Abbruch)
end
%% Animation zeichnen
close all % damit nicht ins falsche Bild gezeichnet wird
if usr_create_animation
  % Zeit-Stützstellen für die Animation vorbereiten
  t_Anfahrt = (0:size(Q_t_Anfahrt,1)-1)'/size(Q_t_Anfahrt,1);
  t_plot = t_Anfahrt;
  Q_plot = Q_t_Anfahrt;
  X_plot = X_t_Anfahrt;
  
  maxduration_animation = 20; % Dauer des mp4-Videos in Sekunden (langsam, damit man folgen kann)
  t_Vid = (0:1/20*(t_plot(end)/maxduration_animation):t_plot(end))';
  I_anim = knnsearch( t_plot , t_Vid );
  I_anim = [I_anim; repmat(length(t_plot),15,1)]; % 15 Standbilder (0.5s) mit letztem Wert
  fprintf('Erstelle Animation mit %d/%d Zeit-Stützstellen\n', length(I_anim), length(t_plot));
  % Animation zeichnen
  anim_filename = fullfile(respath, sprintf('%s_IK_Assembly_Opt_%s',Name, ikobj_all{i_ikobj}));
  s_anim = struct( 'mp4_name', [anim_filename,'.mp4']);%, 'resolution', 200); % 200dpi gibt ungefähr 1080p Breite nach Zuschnitt.
  % Keine Koordinatensysteme plotten. Macht das Bild zu unübersichtlich.
  s_plot = struct( 'ks_legs', [], ...
    'ks_platform', R.NLEG+2, 'straight', 1, 'mode', 4);
  figure(1);clf;hold all;
  set(1, 'name', sprintf('Anim'), ...
    'color','w', 'NumberTitle', 'off', 'units','normalized',...
    'outerposition',[0 0 1 1]); % Vollbild, damit GIF größer wird
  view(3); axis auto; hold on;
  % Keine Achsbeschriftungen
  grid off;
  view(3);
  set(gca,'XTICKLABEL',{}, 'YTICKLABEL', {}, 'ZTICKLABEL',{});
  set(gca,'xtick',[],'ytick',[],'ztick',[])
  set(get(gca, 'XAxis'), 'visible', 'off');
  set(get(gca, 'YAxis'), 'visible', 'off');
  set(get(gca, 'ZAxis'), 'visible', 'off');
  set(gca, 'Box', 'off');
  R.anim( Q_plot(I_anim,:), X_plot(I_anim,:), s_anim, s_plot);
end
end