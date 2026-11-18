function plot_map_climatology(parameter,sit_folder,outfolder,version)

%% Plot maps of parameter and anomalies over one annual cycle

%%% Parameter %%%
% 'radar_freeboard'
% 'sea_ice_freeboard'
% 'sea_ice_thickness'
% 'snow_depth'
% 'sea_ice_concentration'
% 'sea_ice_type'
% 'sea_ice_density'

% year = 2012;
% parameter = 'sea_ice_thickness';
% sit_folder = fullfile(cd,'year_round_sit_v2p0');
% outfolder = cd;
% version = 'v2p0';

%% Code dependency (add to MATLAB path)
% inpaint_nans
% polarstereo_fwd
% othercolor
% https://se.mathworks.com/matlabcentral/fileexchange/30564-othercolor

%% Satellite data

sit_files = dir(fullfile(sit_folder,['uit_sit_*_' version '.nc']));
sit_names = cellfun(@(grid_x) extract(grid_x,digitsPattern(8)), {sit_files.name}');
tvec_sit = cell2mat(cellfun(@(grid_x) [str2double(grid_x(1:4)) str2double(grid_x(5:6)) str2double(grid_x(7:8))], sit_names, 'uni',0));

% Grid
sit_filename = fullfile(sit_files(1).folder,sit_files(1).name);
grid_def = ncreadatt(sit_filename,'/','grid');
grid_lat = ncread(sit_filename,'latitude');
grid_lon = ncread(sit_filename,'longitude');
[grid_x,grid_y] = polarstereo_fwd(grid_lat,grid_lon,6378137,0.08181919,70,0);

%% Export setup

outdir = ['general_plots_' version];
if isfolder(fullfile(outfolder,outdir))
else
    mkdir(fullfile(outfolder,outdir))
end

%% Prepare satellite data time series

P = NaN(size(grid_lat,1),size(grid_lat,2),length(sit_files));
sic = NaN(size(grid_lat,1),size(grid_lat,2),length(sit_files));
type = NaN(size(grid_lat,1),size(grid_lat,2),length(sit_files));
clear t_lim
for i = 1:length(sit_files)
    
    % fprintf(['Date = ' num2str(tvec_sit(i,:)) ' \n\n'])

    %% Load sit grids
    
    sit_filename = fullfile(sit_files(i).folder,sit_files(i).name);

    sit_t = datetime(ncread(sit_filename,'time'), 'ConvertFrom', 'epochtime');
    sit_t_interval = ncread(sit_filename,'time_interval');
    t_lim(i,:) = [sit_t - days(floor(sit_t_interval/2)) sit_t + days(floor(sit_t_interval/2))];
    
    P(:,:,i) = ncread(sit_filename,parameter);
    
    sic(:,:,i) = ncread(sit_filename,'sea_ice_concentration');
    type(:,:,i) = ncread(sit_filename,'sea_ice_type');
    
end

if strcmp(parameter,'sea_ice_concentration')
    P = sic;
elseif strcmp(parameter,'sea_ice_type')
    P = type;
end


%% Climatologies

tvec_clim = uniquetol(tvec_sit(:,2:3),1/max(tvec_sit(:,3)),'highest','Byrows',true);

P_clim = NaN(size(grid_x,1),size(grid_x,2),size(tvec_clim,1));
sic_clim = NaN(size(grid_x,1),size(grid_x,2),size(tvec_clim,1));
type_clim = NaN(size(grid_x,1),size(grid_x,2),size(tvec_clim,1));
for i = 1:size(tvec_clim,1)
    N = sum(~isnan(P(:,:,ismembertol(tvec_sit(:,2:3),tvec_clim(i,:),1/max(tvec_sit(:,3)),'Byrows',true))),3);
    temp = nanmean(P(:,:,ismembertol(tvec_sit(:,2:3),tvec_clim(i,:),1/max(tvec_sit(:,3)),'Byrows',true)),3);
    temp(N<5) = NaN;
    P_clim(:,:,i) = temp;
    
    N = sum(~isnan(sic(:,:,ismembertol(tvec_sit(:,2:3),tvec_clim(i,:),1/max(tvec_sit(:,3)),'Byrows',true))),3);
    temp = nanmean(sic(:,:,ismembertol(tvec_sit(:,2:3),tvec_clim(i,:),1/max(tvec_sit(:,3)),'Byrows',true)),3);
    temp(N<5) = NaN;
    sic_clim(:,:,i) = temp;
    
    N = sum(~isnan(type(:,:,ismembertol(tvec_sit(:,2:3),tvec_clim(i,:),1/max(tvec_sit(:,3)),'Byrows',true))),3);
    temp = nanmean(type(:,:,ismembertol(tvec_sit(:,2:3),tvec_clim(i,:),1/max(tvec_sit(:,3)),'Byrows',true)),3);
    temp(N<5) = NaN;
    type_clim(:,:,i) = temp;
end

% % Numbers:
% % climatological mean thickness
% figure; plot(datetime([ones(size(tvec_clim,1),1)*2020 tvec_clim]),permute(nanmean(P_clim,[1 2]),[3 1 2]))
% % interannual std for date 
% std(permute(nanmean(P(:,:,tvec_sit(:,2)==5 & tvec_sit(:,3)==23),[1 2]),[3 1 2])) % may 23
% std(permute(nanmean(P(:,:,tvec_sit(:,2)==9 & tvec_sit(:,3)==7),[1 2]),[3 1 2])) % sep 07



%% Anomalies

P_anom = NaN(size(P));
sic_anom = NaN(size(P));
for i = 1:size(P,3)
    P_anom(:,:,i) = P(:,:,i) - P_clim(:,:,ismembertol(tvec_clim,tvec_sit(i,2:3),1/max(tvec_sit(:,3)),'Byrows',true));
    sic_anom(:,:,i) = sic(:,:,i) - sic_clim(:,:,ismembertol(tvec_clim,tvec_sit(i,2:3),1/max(tvec_sit(:,3)),'Byrows',true));
    type_anom(:,:,i) = type(:,:,i) - type_clim(:,:,ismembertol(tvec_clim,tvec_sit(i,2:3),1/max(tvec_sit(:,3)),'Byrows',true));
end


%% Plotting

warning('off','all')

Titles = {'Dec, Jan, Feb','Mar, Apr, May','Jun, Jul, Aug','Sep, Oct, Nov'};

% Intervals for each parameter
if strcmp(parameter,'radar_freeboard')
    c_int = 0:0.01:0.35;
    c_std_int = 0:0.005:0.1;
    cmp1 = flipud(othercolor('Spectral10',64));
    cmp2 = viridis(64);
    Tstr = [upper(extractBefore(parameter,2)) strrep(parameter(2:end),'_',' ') ' [m]'];
    Tstr_std = ' Standard deviation [m]';
elseif strcmp(parameter,'sea_ice_freeboard')
    c_int = 0:0.02:0.5;
    c_std_int = 0:0.001:0.15;
    cmp1 = flipud(othercolor('Spectral10',64));
    cmp2 = viridis(64);      
    Tstr = [upper(extractBefore(parameter,2)) strrep(parameter(2:end),'_',' ') ' [m]'];
    Tstr_std = ' Standard deviation [m]';
elseif strcmp(parameter,'sea_ice_thickness')
    c_int = 0:0.1:4;
    c_std_int = 0:0.025:1;
    cmp1 = flipud(othercolor('Spectral10',64));
    cmp2 = viridis(64);      
    Tstr = [upper(extractBefore(parameter,2)) strrep(parameter(2:end),'_',' ') ' [m]'];
    Tstr_std = ' Standard deviation [m]';
elseif strcmp(parameter,'snow_depth')
    c_int = 0:0.02:0.5;
    c_std_int = 0:0.001:0.15;
    cmp1 = othercolor('Blues9',64);
    cmp2 = viridis(64);      
    Tstr = [upper(extractBefore(parameter,2)) strrep(parameter(2:end),'_',' ') ' [m]'];
    Tstr_std = ' Standard deviation [m]';
elseif strcmp(parameter,'sea_ice_concentration')
    c_int = 0:5:100;
    c_std_int = 0:1:30;
    cmp1 = flipud(othercolor('PuBu9',64));
    cmp2 = viridis(64);      
    Tstr = [upper(extractBefore(parameter,2)) strrep(parameter(2:end),'_',' ') ' [%]'];
    Tstr_std = ' Standard deviation [pp]';
elseif strcmp(parameter,'sea_ice_type')
    c_int = 0:0.05:1;
    c_std_int = 0:0.01:0.3;
    cmp1 = othercolor('Blues9',64);
    cmp2 = viridis(64);      
    Tstr = [upper(extractBefore(parameter,2)) strrep(parameter(2:end),'_',' ')];
    Tstr_std = ' Standard deviation';
elseif strcmp(parameter,'sea_ice_density')
    c_int = 868:2:918;
    c_std_int = 0:0.5:10;
    cmp1 = othercolor('YlGnBu9',64);
    cmp2 = viridis(64);      
    Tstr = [upper(extractBefore(parameter,2)) strrep(parameter(2:end),'_',' ') ' [kg/m^3]'];
    Tstr_std = ' Standard deviation [kg/m^3]';
end


%% Parameter maps

figure(599); clf(599);
set(gcf,'Position',[161        1321         843         415])
tiledlayout(2,4,"TileSpacing","compact","Padding","compact")

idt = [12,1,2; 3,4,5; 6,7,8; 9,10,11];

for i = 1:4
    nexttile
    cla
    ncpolarm('lat',68,'lon',0);
    try
        temp2 = double(nanmean(sic(:,:,ismember(tvec_sit(:,2),idt(i,:))),3)>15);
        temp2(temp2==0 | grid_lat>88) = NaN;
        geoshow(grid_lat,grid_lon,temp2*-0.01,'DisplayType','Surface')
        
        temp = nanmean(P(:,:,ismember(tvec_sit(:,2),idt(i,:))),3);
        temp2 = inpaint_nans(temp,1); temp(grid_lat>88 & isnan(temp)) = temp2(grid_lat>88 & isnan(temp));
        temp(grid_lat>88.5 | temp==0 | nanmean(sic(:,:,ismember(tvec_sit(:,2),idt(i,:))),3)<15) = NaN;
        
        validMask = ~isnan(temp);
        cleanedValidMask = bwareaopen(validMask, 2);
        isolatedPixels = validMask & ~cleanedValidMask;
        temp(isolatedPixels) = NaN;

        contourfm(grid_lat,grid_lon,temp,c_int,'LineStyle','none')
        
       catch
    end
    
    % contourm(grid_lat,grid_lon,type(:,:,ismember(tvec_sit(:,2),idt(i,:))),[0.5 0.5],'linecolor',[0.6 0.6 0.6],'Fill','off','linewidth',1)
    contourm(grid_lat,grid_lon,nanmean(sic(:,:,ismember(tvec_sit(:,2),idt(i,:))),3),[15 15],'k','Fill','off','linewidth',0.5)

    clim([min(c_int)-0.01 max(c_int)])
    colormap(gca,[0.8 0.8 0.8; cmp1])
    textm(61,48.5,Titles{i},'fontsize',12,'fontweight','bold','HorizontalAlignment','right','backgroundcolor','white','margin',0.5)
end
% Add colorbar
clim([min(c_int) max(c_int)])
C = colorbar('eastoutside','fontsize',12,'fontweight','bold','ticks',linspace(min(c_int),max(c_int),5));
C.Label.String = Tstr; C.Label.FontWeight = 'bold'; C.Label.FontSize = 12;
C.Label.Rotation = 270;
pause(0.1);

for i = 1:4
    nexttile
    cla
    ncpolarm('lat',68,'lon',0);
    try
        temp2 = double(nanmean(sic(:,:,ismember(tvec_sit(:,2),idt(i,:))),3)>15);
        temp2(temp2==0 | grid_lat>88) = NaN;
        geoshow(grid_lat,grid_lon,temp2*-0.01,'DisplayType','Surface')
        
        temp = nanstd(P(:,:,ismember(tvec_sit(:,2),idt(i,:))),[],3);
        temp2 = inpaint_nans(temp,1); temp(grid_lat>88 & isnan(temp)) = temp2(grid_lat>88 & isnan(temp));
        temp(grid_lat>88.5 | temp==0 | nanmean(sic(:,:,ismember(tvec_sit(:,2),idt(i,:))),3)<15) = NaN;
        
        validMask = ~isnan(temp);
        cleanedValidMask = bwareaopen(validMask, 2);
        isolatedPixels = validMask & ~cleanedValidMask;
        temp(isolatedPixels) = NaN;

        contourfm(grid_lat,grid_lon,temp,c_std_int,'LineStyle','none')
        
       catch
    end
    
    % contourm(grid_lat,grid_lon,type(:,:,ismember(tvec_sit(:,2),idt(i,:))),[0.5 0.5],'linecolor',[0.6 0.6 0.6],'Fill','off','linewidth',1)
    contourm(grid_lat,grid_lon,nanmean(sic(:,:,ismember(tvec_sit(:,2),idt(i,:))),3),[15 15],'k','Fill','off','linewidth',0.5)

    clim([min(c_std_int)-0.01 max(c_std_int)])
    colormap(gca,[0.8 0.8 0.8; cmp2])
    textm(61,48.5,Titles{i},'fontsize',12,'fontweight','bold','HorizontalAlignment','right','backgroundcolor','white','margin',0.5)
end
% Add colorbar
clim([min(c_std_int) max(c_std_int)])
C = colorbar('eastoutside','fontsize',12,'fontweight','bold','ticks',roundn(linspace(min(c_std_int),max(c_std_int),6),-2));
C.Label.String = Tstr_std; C.Label.FontWeight = 'bold'; C.Label.FontSize = 12;
C.Label.Rotation = 270;
pause(0.1);

% export
set(gcf,'color','w')
outfig1name = [parameter '_seasonal_climatology_' version '.png'];
saveas(gcf,fullfile(outfolder,outdir,outfig1name),'png')


warning('on','all')

end