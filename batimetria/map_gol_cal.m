function map_gol_cal(limites);
% function map_gol_cal(limites);
% Hace el mapa del Golfo de California en la zona limitada por
% limites=[sur norte oeste este]. Las longitudes deben ser negativas. Si los
% limites se omiten se hace un mapa de todo el golfo entre [22 32 -116 -105].
% Hay que modificar las trayectorias para encontrar el archivo con linea de
% costa: costa.mat
% Se utiliza el paquete M Map disponible en http://www2.ocgy.ubc.ca/~rich/
% Modificado del programa de Paula (MLM)

direc=cd;
%Defaults
deflim= [22 32 -116 -105];

if nargin<1
    limites= deflim;  
end    

if isempty(limites)
    limites= deflim;
end


axis([limites(3:4) limites(1:2)])
m_proj('mercator','longitudes',limites(3:4), ...
       'latitudes',limites(1:2));

% costa
%m_gshhs_h('patch',[0.7 0.7 0.7]);
% costa de la regi�n de interes, supuestamente se tarda menos que con
% m_gshhs_h
m_usercoast('bat_surfer_grd','patch',[0.7 0.7 0.7]);
hold on
if(direc(1) == 'c'|direc(1)=='C')
%  tray='c:\manuel\batime\inegi\costa';
else
  tray='d:\manuel\batime\inegi\costa';
end
%load(tray);

%m_plot(xx,yy,'k');


%m_grid('fontsize',14,'lineweight',1,'tickdir','in','box','fancy','linestyle','none');
m_grid('fontsize',14,'lineweight',1,'tickdir','in','box','fancy','linestyle','--');

set(gcf,'renderer','painters')