% bat_serfer.m
% Programa para leer la batimetria que me paso Alberto Amados. Cortarla por
% arriba de 30 grados N ya que parece que no esta bien para latitudes
% mayores (los datos llegan hasta 30.6942). La batimetria tiene una resolucion de 0.0072 en
%  latitud y 0.0062 en longitud (aproximadamente 0.84 km en latitud). Se grafica la batimetria por
% abajo de 30N y se escribe en un archivo *.mat. La batimetria original es
% una matriz de 430 X 500 (latitud X longitud).Los datos positivos son en
% tierra (no se si tengan algun significado los valores positivos).

addpath(genpath('/Users/ameliameggo/Matlab/m_map'));
addpath(genpath('/Users/ameliameggo/Matlab/m_map/'));
savepath;


close all; clear all; dz=1;
arch_out='bat_surfer_grd'; op_out=1;
latmin=28+12/60; latmax=30; lonmin=111.8; lonmax=114.;
limites=[latmin, latmax,-lonmax,-lonmin];



nx=500; ny=430;  % dimensiones batimetria

direc=pwd;
if((direc(1) == 'd')|(direc(1) == 'D'))
  tray(1)='d'; tray_adcp(1)='d'; end


% carga batimetria y corta a la region del mapa.
tray_bat='alby0203_surf_grd.dat';
if(direc(1) == '/'); 
  tray_bat='alby0203_surf_grd.dat';
end
mat_malla=load(tray_bat);
lat_b=mat_malla(:,2); lon_b=mat_malla(:,1); bat=mat_malla(:,3);
lat_b=reshape(lat_b,nx,ny); lon_b=reshape(lon_b,nx,ny); bat=reshape(bat,nx,ny);
lat_b=lat_b'; lon_b=lon_b'; bat=bat';
% elimina valores mayores a 30 grados N
[iy,ix]=find(lat_b <= latmax);
iy1=min(iy); iy2=max(iy); ix1=min(ix); ix2=max(ix);
lat_b=lat_b(iy1:iy2,ix1:ix2); lon_b=lon_b(iy1:iy2,ix1:ix2); bat=bat(iy1:iy2,ix1:ix2);
[ny2,nx2]=size(lat_b);
% grafica mapa
map_gol_cal(limites);
title(['Batimetria con resolucion de 0.0072 X 0.0062 y matriz de ',num2str(ny2,'%3i'),' X ',...
    num2str(nx2,'%3i'),' (lat X lon)'])

  [cs,hlab]=m_contour(lon_b,lat_b,bat,[-200 -400 -600 -800 -1000],'k'); hold on
  n_h=length(hlab);
  for n=1:n_h
    set(hlab(n),'linewidth',1);
  end
  clabel(cs,hlab);
  if(op_out)
    eval(['save ',arch_out,' lon_b lat_b bat'])
  end
%return
