function data = time_series_gate_shp(shp_name,sit_folder,search_radius,outfolder,plot_opt,version)

%% Plot time series of key snow and ice parameters within search radius buffer of gate

% Gate specified as a shapefile line feature

% shp_name = 'Gate_A.shp';
% sit_folder = fullfile(cd,'year_round_sit_v2p0');
% search_radius = 100e3; % in meters
% outfolder = 'test';
% plot_opt = 1; % 0 = no plotting, 1 = plotting
% version = 'v2p0';


%% Code dependency (add to MATLAB path)
% time_series_point
% inpaint_nans
% polarstereo_fwd
% othercolor
% https://se.mathworks.com/matlabcentral/fileexchange/30564-othercolor


%% Load shapefile

S = shaperead(shp_name);
info = shapeinfo(shp_name);

if sum(abs(S.X)>180)>0 && sum(abs(S.Y)>90)>0
    [slat,slon] = projinv(info.CoordinateReferenceSystem,S.X,S.Y);
else
    slat = S.Y;
    slon = S.X;
end

[latb,lonb] = interpm(slat,slon,0.01,'gc');

% Convert to polarstereo
[xb,yb] = polarstereo_fwd(latb,lonb,6378137,0.08181919,70,0);

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
grid_resolution = round(median(diff(grid_x(1,:))));

%% Export setup

outdir = ['general_plots_' version];
if isfolder(fullfile(outfolder,outdir))
else
    mkdir(fullfile(outfolder,outdir))
end

[~,name] = fileparts(shp_name);

%% Identify grid cells within buffer of gate

IDX = rangesearch([grid_x(:) grid_y(:)],[xb yb],search_radius);
IDX = unique([IDX{:}])';

%% Prepare satellite data time series

ifb = NaN(size(grid_lat,1),size(grid_lat,2),length(sit_files));
ifb_unc = NaN(size(grid_lat,1),size(grid_lat,2),length(sit_files));
sd = NaN(size(grid_lat,1),size(grid_lat,2),length(sit_files));
sd_unc = NaN(size(grid_lat,1),size(grid_lat,2),length(sit_files));
sit = NaN(size(grid_lat,1),size(grid_lat,2),length(sit_files));
sit_unc = NaN(size(grid_lat,1),size(grid_lat,2),length(sit_files));
sic = NaN(size(grid_lat,1),size(grid_lat,2),length(sit_files));
type = NaN(size(grid_lat,1),size(grid_lat,2),length(sit_files));
mpf = NaN(size(grid_lat,1),size(grid_lat,2),length(sit_files));
for i = 1:length(sit_files)
    
    % fprintf(['Date = ' num2str(tvec_sit(i,:)) ' \n\n'])

    %% Load sit grids
    
    sit_filename = fullfile(sit_files(i).folder,sit_files(i).name);

    sit(:,:,i) = ncread(sit_filename,'sea_ice_thickness');
    sit_unc(:,:,i) = ncread(sit_filename,'sea_ice_thickness_unc');
    
    ifb(:,:,i) = ncread(sit_filename,'sea_ice_freeboard');
    ifb_unc(:,:,i) = ncread(sit_filename,'sea_ice_freeboard_unc');
    
    sd(:,:,i) = ncread(sit_filename,'snow_depth');
    sd_unc(:,:,i) = ncread(sit_filename,'snow_depth_unc');
    
    sic(:,:,i) = ncread(sit_filename,'sea_ice_concentration');
    type(:,:,i) = ncread(sit_filename,'sea_ice_type');
    % mpf(:,:,i) = ncread(sit_filename,'melt_pond_fraction');
    

end

% Ice draft

sid = sit - ifb;
sid((isnan(sic) | sic<15) & isnan(sid)) = 0;
sid_unc = sit_unc - ifb_unc;
   
    
%% Reshape grids to 2D

ifb = reshape(ifb,[size(grid_lat,1)*size(grid_lat,2) size(tvec_sit,1)]);
ifb_unc = reshape(ifb_unc,[size(grid_lat,1)*size(grid_lat,2) size(tvec_sit,1)]);

sid = reshape(sid,[size(grid_lat,1)*size(grid_lat,2) size(tvec_sit,1)]);
sid_unc = reshape(sid_unc,[size(grid_lat,1)*size(grid_lat,2) size(tvec_sit,1)]);

sd = reshape(sd,[size(grid_lat,1)*size(grid_lat,2) size(tvec_sit,1)]);
sd_unc = reshape(sd_unc,[size(grid_lat,1)*size(grid_lat,2) size(tvec_sit,1)]);

sit = reshape(sit,[size(grid_lat,1)*size(grid_lat,2) size(tvec_sit,1)]);
sit_unc = reshape(sit_unc,[size(grid_lat,1)*size(grid_lat,2) size(tvec_sit,1)]);

sic = reshape(sic,[size(grid_lat,1)*size(grid_lat,2) size(tvec_sit,1)]);
type = reshape(type,[size(grid_lat,1)*size(grid_lat,2) size(tvec_sit,1)]);
% mpf = reshape(mpf,[size(grid_lat,1)*size(grid_lat,2) size(tvec_sit,1)]);

%% Time series at reference location

ifb_t = nanmean(ifb(IDX,:))';
ifb_unc_t = nanmean(ifb_unc(IDX,:))';

sid_t = nanmean(sid(IDX,:))';
sid_unc_t = nanmean(sid_unc(IDX,:))';

sd_t = nanmean(sd(IDX,:))';
sd_unc_t = nanmean(sd_unc(IDX,:))';

sit_t = nanmean(sit(IDX,:))';
sit_unc_t = nanmean(sit_unc(IDX,:))';

sic_t = nanmean(sic(IDX,:))';
type_t = nanmean(type(IDX,:))';
% mpf_t = nanmean(mpf(IDX,:))';

%% Time series per unit distance across gate

% Distance along gate (in direction specified by shapefile coords) in m
D_gate = sqrt((xb - xb(1)).^2 + (yb - yb(1)).^2);

% Grid indices per gate coord
IDX2 = rangesearch([grid_x(:) grid_y(:)],[xb yb],search_radius);

% Time series at gate coords
ifb_gate_t = cell2mat(cellfun(@(x) nanmean(ifb(x,:),1), IDX2, 'uni', 0));
ifb_unc_gate_t = cell2mat(cellfun(@(x) nanmean(ifb_unc(x,:),1), IDX2, 'uni', 0));

sid_gate_t = cell2mat(cellfun(@(x) nanmean(sid(x,:),1), IDX2, 'uni', 0));
sid_unc_gate_t = cell2mat(cellfun(@(x) nanmean(sid_unc(x,:),1), IDX2, 'uni', 0));

sd_gate_t = cell2mat(cellfun(@(x) nanmean(sd(x,:),1), IDX2, 'uni', 0));
sd_unc_gate_t = cell2mat(cellfun(@(x) nanmean(sd_unc(x,:),1), IDX2, 'uni', 0));

sit_gate_t = cell2mat(cellfun(@(x) nanmean(sit(x,:),1), IDX2, 'uni', 0));
sit_unc_gate_t = cell2mat(cellfun(@(x) nanmean(sit_unc(x,:),1), IDX2, 'uni', 0));

sic_gate_t = cell2mat(cellfun(@(x) nanmean(sic(x,:),1), IDX2, 'uni', 0));
type_gate_t = cell2mat(cellfun(@(x) nanmean(type(x,:),1), IDX2, 'uni', 0));
% mpf_gate_t = cell2mat(cellfun(@(x) nanmean(mpf(x,:),1), IDX2, 'uni', 0));


%% Climatologies

tvec_clim = uniquetol(tvec_sit(:,2:3),1/max(tvec_sit(:,3)),'highest','Byrows',true);
[~,LocB] = ismembertol(tvec_sit(:,2:3),tvec_clim,1/max(tvec_sit(:,3)),'Byrows',true);

ifb_t_clim = accumarray(LocB,ifb_t(:,1),size(tvec_clim(:,1)),@nanmean);
ifb_unc_t_clim = accumarray(LocB,ifb_unc_t(:,1),size(tvec_clim(:,1)),@nanmean);

sid_t_clim = accumarray(LocB,sid_t(:,1),size(tvec_clim(:,1)),@nanmean);
sid_unc_t_clim = accumarray(LocB,sid_unc_t(:,1),size(tvec_clim(:,1)),@nanmean);

sd_t_clim = accumarray(LocB,sd_t(:,1),size(tvec_clim(:,1)),@nanmean);
sd_unc_t_clim = accumarray(LocB,sd_unc_t(:,1),size(tvec_clim(:,1)),@nanmean);

sit_t_clim = accumarray(LocB,sit_t(:,1),size(tvec_clim(:,1)),@nanmean);
sit_unc_t_clim = accumarray(LocB,sit_unc_t(:,1),size(tvec_clim(:,1)),@nanmean);

sic_t_clim = accumarray(LocB,sic_t(:,1),size(tvec_clim(:,1)),@nanmean);
type_t_clim = accumarray(LocB,type_t(:,1),size(tvec_clim(:,1)),@nanmean);


%% Calculate anomalies

ifb_t_anom = ifb_t(:,1) - ifb_t_clim(LocB);
sid_t_anom = sid_t(:,1) - sid_t_clim(LocB);
sd_t_anom = sd_t(:,1) - sd_t_clim(LocB);
sit_t_anom = sit_t(:,1) - sit_t_clim(LocB);
sic_t_anom = sic_t(:,1) - sic_t_clim(LocB);
type_t_anom = type_t(:,1) - type_t_clim(LocB);

%% Final data

data = [];
data.tvec = tvec_sit;
data.ifb_t = ifb_t;
data.ifb_unc_t = ifb_unc_t;
data.ifb_t_anom = ifb_t_anom;
data.sid_t = sid_t;
data.sid_unc_t = sid_unc_t;
data.sid_t_anom = sid_t_anom;
data.sd_t = sd_t;
data.sd_unc_t = sd_unc_t;
data.sd_t_anom = sd_t_anom;
data.sit_t = sit_t;
data.sit_unc_t = sit_unc_t;
data.sit_t_anom = sit_t_anom;
data.sic_t = sic_t;
data.sic_t_anom = sic_t_anom;
data.type_t = type_t;
data.type_t_anom = type_t_anom;

% For all individual grid cells along gate
data.lat_idx = grid_lat(IDX)';
data.lon_idx = grid_lon(IDX)';
data.ifb_t_idx = ifb(IDX,:)';
data.ifb_unc_t_idx = ifb_unc(IDX,:)';
data.sid_t_idx = sid(IDX,:)';
data.sid_unc_t_idx = sid_unc(IDX,:)';
data.sd_t_idx = sd(IDX,:)';
data.sd_unc_t_idx = sd_unc(IDX,:)';
data.sit_t_idx = sit(IDX,:)';
data.sit_unc_t_idx = sit_unc(IDX,:)';
data.sic_t_idx = sic(IDX,:)';
data.type_t_idx = type(IDX,:)';


if plot_opt==0
else

    %% Plotting
    
    idt = zeros(size(tvec_sit,1),1); idt(tvec_sit(:,2)>4 & tvec_sit(:,2)<10) = 1;
    
    idn = ~isnan(ifb_t) | ~isnan(sit_t);
    
    %% Time series
    
    figure(99); clf(99);
    set(gcf,'Position',[1	41	1280	1317.60000000000])
    tiledlayout(4,1,"TileSpacing","compact","Padding","compact")
    
    cmap1 = othercolor('Blues5',5);
    nexttile
    cla
    hold on
    patch([datenum(tvec_sit(idn,:)); flipud(datenum(tvec_sit(idn,:)))],[ifb_t(idn); flipud(-sid_t(idn))],cmap1(5,:),'facealpha',0.5,'edgecolor','none')
    patch([datenum(tvec_sit(idn,:)); flipud(datenum(tvec_sit(idn,:)))],[ifb_t(idn) + sd_t(idn); flipud(ifb_t(idn))],cmap1(3,:),'facealpha',0.5,'edgecolor','none')
    grid on
    box on
    datetick('x','yyyy')
    xlim([datenum(tvec_sit(1,:)) datenum(tvec_sit(end,:))])
    % xlabel('Year')
    xticklabels('')
    ylim([-4 2])
    title(strrep(name,'_',' '))
    ylabel('Relative Height ASL [m]')
    legend('Sea Ice Thickness','Snow Thickness','Location','northeast')
    set(gca,'fontsize',14)
    hold off
    
    cmap1 = othercolor('PRGn4',2);
    nexttile
    cla
    hold on
    plot(datenum(tvec_sit(idn,:)),sd_t_anom(idn),'color',cmap1(1,:),'linewidth',2)
    plot(datenum(tvec_sit(idn,:)),ifb_t_anom(idn),'color',cmap1(2,:),'linewidth',2)
    grid on
    box on
    datetick('x','yyyy')
    xlim([datenum(tvec_sit(1,:)) datenum(tvec_sit(end,:))])
    % xlabel('Year')
    xticklabels('')
    ylim([-0.3 0.3])
    ylabel('Thickness Anomaly [m]')
    legend('Snow Depth Anomaly','Sea Ice Freeboard Anomaly','Location','northeast')
    set(gca,'fontsize',14)
    hold off
    
    cmap1 = othercolor('RdBu4',2);
    nexttile
    cla
    hold on
    plot(datenum(tvec_sit(idn,:)),sit_t_anom(idn),'color',cmap1(1,:),'linewidth',2)
    grid on
    box on
    datetick('x','yyyy')
    xlim([datenum(tvec_sit(1,:)) datenum(tvec_sit(end,:))])
    % xlabel('Year')
    xticklabels('')
    ylim([-2 2])
    ylabel('Thickness Anomaly [m]')
    legend('Sea Ice Thickness Anomaly','Location','northeast')
    set(gca,'fontsize',14)
    hold off
    
    cmap1 = othercolor('BuOr_12',5);
    nexttile
    cla
    hold on
    plot(datenum(tvec_sit(idn,:)),sic_t(idn),'color',cmap1(1,:),'linewidth',2)
    plot(datenum(tvec_sit(idn,:)),type_t(idn)*100,'color',cmap1(4,:),'linewidth',2)
    grid on
    box on
    datetick('x','yyyy')
    xlim([datenum(tvec_sit(1,:)) datenum(tvec_sit(end,:))])
    xlabel('Year')
    ylim([0 100])
    ylabel('Percentage [%]')
    legend('Sea Ice Concentration','Fraction of MYI','Location','southeast')
    set(gca,'fontsize',14)
    hold off
    
    ax1 = nexttile(1);
    ax2 = nexttile(2);
    ax3 = nexttile(3);
    ax4 = nexttile(4);
    linkaxes([ax1 ax2 ax3 ax4],'x')

    cmap1 = othercolor('Blues5',5);
    ax5 = axes('Position',[.07 .07 .1 .1]);
    ncpolarm('lat',floor(min(grid_lat(IDX))),'lon',-45)
    scatterm(grid_lat(IDX),grid_lon(IDX),5,cmap1(3,:),'filled')
    plotm(latb,lonb,'color',cmap1(5,:),'linewidth',2)
    
    % export
    set(gcf,'color','w')
    outfig1name = ['time_series_' name '_' version '.png'];
    saveas(gcf,fullfile(outfolder,outdir,outfig1name),'png')


    %% Time series across gate

    % plot wrt lat or lon?
    if range(latb) > 3
        spanb = latb;
        yaxname = 'Latitude';
    else
        spanb = lonb;
        yaxname = 'Longitude';
    end

    figure(100); clf(100);
    set(gcf,'Position',[1	41	1280	1317.60000000000])
    tiledlayout(4,1,"TileSpacing","compact","Padding","compact")
    
    cmap1 = flipud(othercolor('PuBu9',64));
    nexttile
    cla
    hold on
    imagesc(datenum(tvec_sit),spanb(1:end-1),ifb_gate_t(1:end-1,:))
    box on
    set(gca,'layer','top')
    datetick('x','yyyy')
    xlim([datenum(tvec_sit(1,:)) datenum(tvec_sit(end,:))])
    % xlabel('Year')
    xticklabels('')
    ylim([min(spanb) max(spanb)])
    title(strrep(name,'_',' '))
    ylabel(yaxname)
    clim([0 0.5])
    colormap(gca,flipud(cmap1))
    c = colorbar;
    c.Label.String = 'Sea Ice Freeboard [m]';
    set(gca,'fontsize',14)
    hold off
    
    nexttile
    cla
    hold on
    imagesc(datenum(tvec_sit),spanb(1:end-1),sd_gate_t(1:end-1,:))
    box on
    set(gca,'layer','top')
    datetick('x','yyyy')
    xlim([datenum(tvec_sit(1,:)) datenum(tvec_sit(end,:))])
    % xlabel('Year')
    xticklabels('')
    ylim([min(spanb) max(spanb)])
    ylabel(yaxname)
    clim([0 0.5])
    colormap(gca,flipud(cmap1))
    c = colorbar;
    c.Label.String = 'Snow Depth [m]';
    set(gca,'fontsize',14)
    hold off
    
    cmap2 = flipud(othercolor('YlGnBu9',64));
    nexttile
    cla
    hold on
    imagesc(datenum(tvec_sit),spanb(1:end-1),sit_gate_t(1:end-1,:))
    box on
    set(gca,'layer','top')
    datetick('x','yyyy')
    xlim([datenum(tvec_sit(1,:)) datenum(tvec_sit(end,:))])
    % xlabel('Year')
    xticklabels('')
    ylim([min(spanb) max(spanb)])
    ylabel(yaxname)
    clim([0 4])
    colormap(gca,flipud(cmap2))
    c = colorbar;
    c.Label.String = 'Sea Ice Thickness [m]';
    set(gca,'fontsize',14)
    hold off
    
    cmap3 = flipud(othercolor('BrBG9',128));
    cmap3 = cmap3(65:end,:);
    nexttile
    cla
    hold on
    imagesc(datenum(tvec_sit),spanb(1:end-1),type_gate_t(1:end-1,:)*100)
    box on
    set(gca,'layer','top')
    datetick('x','yyyy')
    xlim([datenum(tvec_sit(1,:)) datenum(tvec_sit(end,:))])
    xlabel('Year')
    % xticklabels('')
    ylim([min(spanb) max(spanb)])
    ylabel(yaxname)
    clim([0 100])
    colormap(gca,cmap3)
    c = colorbar;
    c.Label.String = 'Fraction of MYI [%]';
    set(gca,'fontsize',14)
    hold off
    
    ax1 = nexttile(1);
    ax2 = nexttile(2);
    ax3 = nexttile(3);
    ax4 = nexttile(4);
    linkaxes([ax1 ax2 ax3 ax4],'x')

    cmap1 = othercolor('Blues5',5);
    ax5 = axes('Position',[.07 .07 .1 .1]);
    ncpolarm('lat',floor(min(grid_lat(IDX))),'lon',-45)
    % scatterm(grid_lat(IDX),grid_lon(IDX),5,cmap1(3,:),'filled')
    plotm(latb,lonb,'color',cmap1(5,:),'linewidth',2)
    
    % Update font size
    allAxes = findobj(gcf, 'type', 'axes');
    set(allAxes, 'FontSize', 14);
    
    % export
    set(gcf,'color','w')
    outfig2name = ['time_series_across_gate_' name '_' version '.png'];
    saveas(gcf,fullfile(outfolder,outdir,outfig2name),'png')

end

end
