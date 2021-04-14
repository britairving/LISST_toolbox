function [cfg, data_proc,meta_proc] = LISST_data_qaqc_expert(cfg,data_proc,meta_proc)
%function LISST_data_qaqc_expert runs quality checks on data
%
%  Syntax:
%    [cfg, data_proc] = LISST_data_qaqc_auto(cfg,data_proc)
%
%  Description:
%    Performs automated quality control on processed LISST data.
%    Performs qc tests descripted in manual, as well as test included in
%    InLineAnaylsis processLISST.m script.
%
%  Refereces:
%    LISST-Deep-Users-Manual-May-2013.pdf "STEP BY STEP PROCEDURE: DATA QUALITY CONTROL"
%    https://seabass.gsfc.nasa.gov/archive/MAINE/boss/EXPORTS/exportsnp/documents/EXPORTS-EXPORTSNP_InLine-LISST-Processing_R1.pdf
%
%  Notes:
%    Incorporate in the future...
%      http://www.sequoiasci.com/article/how-accurate-is-the-concentration-measurement-by-lisst-instruments/
%      https://www.sequoiasci.com/article/the-influence-of-particles-outside-the-size-range-of-the-lisst/
%      https://www.sequoiasci.com/article/interpreting-particle-data/
%    Review and incorporate ambient light effects Andrews et al 2011 https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2010WR009841
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>

% IV. THE REASONABLENESS TEST
%
% Ever hear an old guy in an audience say, with great authority: “I don’t
% believe it!”. The final test of your results is reasonableness. This test
% uses all the knowledge and wisdom of your past experience to evaluate the
% current result. If it does not make sense, there is probably a good
% reason.
%
% Here are some criteria that you can use to evaluate if ‘it’ does make sense:
%
%     Is the net scattering smooth? The final net scattering variable cscat
%     should be smooth. The nature of light scattering by particles is such
%     that sudden spikes in an individual ring can not arise. For example,
%     a cscat curve for the 123rd scan shows a spike at ring 23. This is
%     not possible. No particle size produces a spike at a single ring.
%     What is the likely cause? Probably the background data is not good,
%     or the estimate of optical transmission is not good. Sometimes, a
%     minor adjustment of the background or the transmission can remove
%     such a spike. At other times, you can use the criterion of a smooth
%     cscat and sophisticated mathematical routines to find a best fit
%     estimate of transmission. We have successfully used a minimization
%     technique that seeks the smoothest cscat by varying transmission.
%     Such routines are available in Matlab. The function fmins is one such
%     function.
%     A persistent high value at the inner rings spells trouble. In
%     situations where large particles are not always present, a persistent
%     high value at the inner rings can not exist in cscat. Again, this
%     would point to a bad estimate of the background zscat.

close all
%% 0 | Set up script basics
fprintf('NEED TO SET THIS UP!!!\n')
keyboard
all_casts     = unique(data_proc.cast)';
flagged_scans = zeros(size(data_proc.date));
flag_value    = NaN; % NaN or 0 - does change how automatic outlier tests perform


% don't overwrite already flagged values in case different flag codes
idx_alreadyflagged = data_proc.quality_flag(isuspect) ~= flag.not_evaluated;
isuspect = isuspect(idx_alreadyflagged);
data_proc.quality_flag(isuspect) = flag.questionable;

idx = data_proc.quality_flag ~= 2;
%% Loop through each cast and identify outliers
makefig; 
ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);

for ncast = all_casts
  %% pull out subset with cast data
  cla(ax1); cla(ax2);
  %% pull out subset with cast data
  idx_cast = find(data_proc.cast == ncast);
  cscat = data_proc.cscat(idx_cast,:);
  VSD   = data_proc.VSD(idx_cast,:);
  
  %% Make sure proper flag value - will change automated qc performance..
  if isnan(flag_value)
    cscat(cscat == 0) = flag_value;
  elseif flag_value == 0
    cscat(isnan(cscat)) = flag_value;
  end
  
  %% remove data already flagged
  idx_flagged = data_proc.quality_flag(idx_cast) ~= flag.good & data_proc.quality_flag(idx_cast) ~= flag.not_evaluated;
  cscat(idx_flagged,:) = flag_value;
  VSD(idx_flagged,:)   = flag_value;
 
  %% Test | profiles with tail
  idx_flagged = VSD < 10e-70;
  idx_flagged = any(idx_flagged,2);
  %[VSD,cscat] = highlight_and_choose_to_remove(idx_flagged,VSD,cscat,'PSD<10e-100');
  cscat(idx_flagged,:) = flag_value;
  VSD(idx_flagged,:)   = flag_value;
  data_proc.quality_flag(idx_cast(idx_flagged)) = flag.bad;
  flagged_scans(idx_cast(idx_flagged)) = 1;
  
  %% Move to next cast if all values are zero (or NaN)
  if all(all(isnan(cscat))) || all(all(cscat == 0)) || sum(sum(isfinite(cscat))) < 100
    continue
  end

  %% plot uncorrected scatter
  h1 = plot(ax1,1:32,cscat,'-');
  hold(ax1,'on'); grid(ax1,'on')
  try
    lims = prctile(cscat(:,10:end),[1 98]);
    ax1.YLim = [0 max(max(lims))]; %[min(min(lims)) max(max(lims))]; % 1st and 99th percentiles
  end
  ax1.YLabel.String = 'Corrected Scattering';
  ax1.XLabel.String = 'Ring Detector No.';
  ax1.XLim = [1 32]; ax1.XTick = [1:32];
  title(ax1,{[data_proc.datfile{idx_cast(1)} ' | Cast ' num2str(ncast)]; ['Background: ' strrep(data_proc.zscfile{idx_cast(1)},'_','\_')]});
  
  h2 = plot(ax2,cfg.dias,VSD,'-');
  hold(ax2,'on'); grid(ax2,'on')
  ax2.XLim = [floor(min(cfg.dias)) ceil(max(cfg.dias)+10)];
  ax2.YLim(2) = 1; ax2.XScale = 'log'; ax2.YScale = 'log';
  ax2.YLabel.String = 'Calibrated VSD [\muL/L]';
  ax2.XLabel.String = 'Mean Diameter [\mum]';
  ax2.XTick = cfg.dias; ax2.XAxis.TickLabelFormat = '%.1f';
  ax2.XAxis.TickLabelRotation = 45;
  ax2.XLim = [cfg.dias(1) cfg.dias(end)]; ax2.XDir = 'reverse'; % ring #1 == largest size bin
  
  done = 0;
  while ~done
    if userprompted_qc
      fprintf('\n')
      fprintf('How does the data from cast %d look?\n',ncast)
      fprintf('   <1>  All questionable  | flag whole cast\n')
      fprintf('   <2>  Some questionable | select questionable data\n')
      fprintf('   <3>  All good          | do not flag anything\n')
      fprintf('   <99> Stop              | keyboard mode\n')
      chc = input('   Enter choice: ');
    else
      chc = 2;
    end
    % default to all good
    if isempty(chc)
      chc = 3; 
    end 
    switch chc
      case 1
        done = 1;
        % don't overwrite the previously flagged data
        idx_flagtherest = ~ismember(idx_cast,flagged_scans);
        data_proc.quality_flag(idx_cast(idx_flagtherest)) = flag.questionable;
        flagged_scans(idx_cast(idx_flagtherest)) = 1;
      case 2
        %% Test | Remove scans with multiple peaks (spiky scans)
        idx_flagged = [];
        n_peaks     = [];
        for nscan = 1:size(VSD,1)
          if all(isnan(cscat(nscan,:)))
            continue
          end
          [pks,~,~,~] = findpeaks(cscat(nscan,:),'threshold',3);
          if numel(pks) > 2 %~isempty(pks) %  numel(pks) > 1
            idx_flagged = [idx_flagged; nscan];
            n_peaks = [n_peaks;  numel(pks)];
          end
          %fprintf('Depth = %.1fm | Scan # %d\n',data_proc.depth(idx_cast(nscan)),nscan)
        end
        if userprompted_qc
          [VSD,cscat] = highlight_and_choose_to_remove(idx_flagged,VSD,cscat,'Smoothness test');
        else
          % return if no flagged data
          if all(idx_flagged == 0 | idx_flagged == 1) && ~islogical(idx_flagged)
            idx_flagged = logical(idx_flagged);
          end
          fprintf('%s\n','Smoothness test')
          if sum(idx_flagged) > 0
            hflag1 = plot(ax1,1:32,cscat(idx_flagged,:),'r-','LineWidth',2);
            hflag2 = plot(ax2,cfg.dias,VSD(idx_flagged,:),'r-','LineWidth',2);
            htext  = text(ax1,0.1,0.85,['Data flagged using test: Smoothness test'],'units','normalized','Color','r','fontweight','bold');
            %figname = fullfile(cfg.path.figs,[cfg.project '_qc_Cast' num2str(ncast)]);
            %standard_printfig_lowrespng(figname)
          end

          % don't overwrite the previously flagged data
          idx_flagtherest = ~ismember(idx_cast,flagged_scans);
          data_proc.quality_flag(idx_cast(idx_flagtherest)) = flag.questionable;
          flagged_scans(idx_cast(idx_flagtherest)) = 1;
        end
        %         %% Test | Iterative median outliers test
        %         % First - remove whole scan
        %         idx_flagged = isoutlier(cscat,'median','ThresholdFactor',7); % true if more than ThresholdFactor scaled MAD from median. The scaled MAD is defined as c*median(abs(A-median(A))), where c=-1/(sqrt(2)*erfcinv(3/2)).
        %         idx_flagged = all(idx_flagged,2);      %% ALL
        %         if userprompted_qc
        %           [VSD,cscat] = highlight_and_choose_to_remove(idx_flagged,VSD,cscat,'whole scan median outlier');
        %         else
        %           % don't overwrite the previously flagged data
        %           idx_flagtherest = ~ismember(idx_cast,flagged_scans);
        %           data_proc.quality_flag(idx_cast(idx_flagtherest)) = flag.questionable;
        %           flagged_scans(idx_cast(idx_flagtherest)) = 1;
        %         end
        %         % now loop through and remove outliers
        %         for npass = 1:2
        %           idx_flagged = isoutlier(cscat,'median','ThresholdFactor',33); % true if more than ThresholdFactor scaled MAD from median. The scaled MAD is defined as c*median(abs(A-median(A))), where c=-1/(sqrt(2)*erfcinv(3/2)).
        %           %idx_flagged = all(idx_flagged,2);      %% ALL
        %           %idx_flagged = any(idx_flagged,2);      %% ANY
        %           idx_flagged = sum(idx_flagged,2) > 10;   %% SOME
        %           if userprompted_qc
        %             [VSD,cscat] = highlight_and_choose_to_remove(idx_flagged,VSD,cscat,['test #' num2str(npass) ' median outliers']);
        %           else
        %             % don't overwrite the previously flagged data
        %             idx_flagtherest = ~ismember(idx_cast,flagged_scans);
        %             data_proc.quality_flag(idx_cast(idx_flagtherest)) = flag.questionable;
        %             flagged_scans(idx_cast(idx_flagtherest)) = 1;
        %           end
        %         end
        done = 1;
      case 3
        done = 1;
      case 99
        fprintf('...in keyboard mode, enter "dbcont" to continue\n')
        keyboard
        done = 0;
      otherwise
        fprintf('...unknown selection, try again\n')
        done = 0;
    end
  end % Choose what to do with cast
end

%% 5 | Save qc'd data to .mat file
% if ~ cfg.testing % only save if not in testing mode
fprintf('Saving LISST QC data to %s\n',cfg.path.file_qc)
if exist(cfg.path.file_qc,'file')
  movefile(cfg.path.file_qc, cfg.path.oldorunused)
end
save(cfg.path.file_qc,'cfg','data_proc');
% end

%% FUNCTION highlight_and_choose_to_remove
  function [VSD,cscat] = highlight_and_choose_to_remove(idx_flagged,VSD,cscat,test_string)
    % return if no flagged data
    if all(idx_flagged == 0 | idx_flagged == 1) && ~islogical(idx_flagged)
      idx_flagged = logical(idx_flagged);
    end
    fprintf('%s\n',test_string)
    if sum(idx_flagged) == 0
      return
    end
    hflag1 = plot(ax1,1:32,cscat(idx_flagged,:),'r-','LineWidth',2);
    hflag2 = plot(ax2,cfg.dias,VSD(idx_flagged,:),'r-','LineWidth',2);
    htext  = text(ax1,0.1,0.85,['Data flagged using test: ' test_string],'units','normalized','Color','r','fontweight','bold');
    
    fprintf('Do you want to remove the flagged data? Default is yes (1)\n')
    fprintf('  <0>  No        \n')
    fprintf('  <1>  Yes (default)\n')
    fprintf('  <99> Stop        \n')
    rm_chc = input(' Enter choice: ');
    if isempty(rm_chc); rm_chc = 1; end
    if rm_chc == 1
      cscat(idx_flagged,:) = flag_value;
      VSD(idx_flagged,:)   = flag_value;
      data_proc.quality_flag(idx_cast(idx_flagged)) = flag.bad;
      flagged_scans(idx_cast(idx_flagged)) = 1;
    end
    if rm_chc == 99
      keyboard
    end
    delete(hflag1);
    delete(hflag2);
    delete(h1);
    delete(h2);
    delete(htext);
    h1 = plot(ax1,1:32,cscat,'-');
    h2 = plot(ax2,cfg.dias,VSD,'-');
    %ax1.YLim = 'auto';
  end %% FUNCTION highlight_and_choose_to_remove

%% function plot_flagged_data
%   function plot_flagged_data(data_proc,idx_cast,cfg)
%     clr_questionable         = 'y';%color to plot questionable data
%     axcolor                  = [0.8 0.8 0.8];
%     dntick = 'HH:MM';
%     dc = data_proc(idx_cast);
%     % pull out DAT and zscat file
%     [~,datfile1,ext] = fileparts(dc.datfile{1});
%     datfile = strrep([datfile1 ext],'_','\_');
%     [~,zscfile,ext]  = fileparts(dc.zscfile{1});
%     zscfile = strrep([zscfile ext],'_','\_');
%
%     % pull out where data is flagged as bad and questionable
%     ibad = dc.quality_flag == flag.bad;
%     ique = dc.quality_flag == flag.questionable;
%
%     % convert transmission to percent
%     t = dc.tau*100;
%     % make figure
%     fig = makefig; fig.Name = datfile;
%     ax1 = subplot(2,3,1);
%     ax2 = subplot(2,3,4);
%     ax3 = subplot(2,3,[2 5]);
%     ax4 = subplot(2,3,[3 6]);
%
%     %% plot vs time
%     % Laser reference value
%     plot(ax1,dc.date, dc.laserPowerEnteringWater,'k*','DisplayName','Laser Reference');
%     ax1.YLabel.String = 'laser reference [mW]'; hold(ax1,'on'); grid(ax1,'on'); datetick(ax1,'x',dntick,'keeplimits','keepticks');
%     plot(ax1,dc.date(ibad), dc.laserPowerEnteringWater(ibad),'r*','DisplayName','bad','MarkerSize',4)
%     plot(ax1,dc.date(ique), dc.laserPowerEnteringWater(ique),'*','color',clr_questionable,'DisplayName','questionable','MarkerSize',4)
%     %legend(ax1,'show','location','se')
%     title(ax1,{['data: ' datfile]; ['zscat: ' zscfile]})
%     % Transmission
%     plot(ax2,dc.date, t,'k*','DisplayName','transmission');
%     ax2.YLabel.String = 'transmission [%]'; hold(ax2,'on'); grid(ax2,'on'); datetick(ax2,'x',dntick,'keeplimits','keepticks');
%     ax2.YLim = [0 110];
%     plot(ax2,dc.date(ibad), t(ibad),'r*','DisplayName','bad','MarkerSize',4)
%     plot(ax2,dc.date(ique), t(ique),'*','color',clr_questionable,'DisplayName','questionable','MarkerSize',4);
%     ax2.XLabel.String = 'Date';
%     %legend(ax2,'show','location','se')
%     %% plot vs depth
%     c = zeros(size(t,1),3);% black [0 0 0]
%     c(ibad,1) = 1; % red [1 0 0]
%     if strcmp(clr_questionable,'y')
%       c(ique,1:2) = 1; % yellow [1 1 0]
%     else
%       c(ique,3)   = 1; % blue [0 0 1]
%     end
%     scatter(ax3,dc.date,dc.depth,35,c,'filled')
%     ax3.YDir = 'reverse'; datetick(ax3,'x',dntick,'keeplimits','keepticks'); ax3.XLabel.String = 'Date';
%     ax3.YLabel.String = 'Depth [m]'; grid(ax3,'on')
%     text(ax3,0.05,0.20,'good','Color','k','units','normalized','fontweight','bold','HorizontalAlignment','left')
%     text(ax3,0.05,0.15,'bad','Color','r','units','normalized','fontweight','bold','HorizontalAlignment','left')
%     text(ax3,0.05,0.10,'questionable','Color',clr_questionable,'units','normalized','fontweight','bold','HorizontalAlignment','left')
%     title(ax3,[datestr(dc.date(1),'dd-mmm-yyyy HH:MM') ' to ' datestr(dc.date(end),'HH:MM')])
%     % transmission vs depth
%     scatter(ax4,t,dc.depth,35,c,'filled')
%     ax4.YDir = 'reverse';
%     ax4.YLabel.String = 'Depth [m]'; grid(ax4,'on'); ax4.XLabel.String = 'Transmission [%]';
%     text(ax4,0.05,0.20,'good','Color','k','units','normalized','fontweight','bold','HorizontalAlignment','left')
%     text(ax4,0.05,0.15,'bad','Color','r','units','normalized','fontweight','bold','HorizontalAlignment','left')
%     text(ax4,0.05,0.10,'questionable','Color',clr_questionable,'units','normalized','fontweight','bold','HorizontalAlignment','left')
%     title(ax4,[datestr(dc.date(1),'dd-mmm-yyyy HH:MM') ' to ' datestr(dc.date(end),'HH:MM')])
%     % change axe color so easier to see
%     ax1.Color = axcolor;
%     ax2.Color = axcolor;
%     ax3.Color = axcolor;
%     ax4.Color = axcolor;
%     figname = fullfile(cfg.path.figs,[cfg.project '_' datfile1 '_qc' cfg.zscat_choice]);
%     standard_printfig_lowrespng(figname)
%     close
%    
%   end %% function plot_flagged_data

end %% MAIN FUNCTION

%% POSSIBLE CAUSES
% <http://www.sequoiasci.com/article/interpreting-particle-data/>
% OPTICAL TRANSMISSION: Check the optical transmission. If below 30% the
% particle size distribution (PSD) may show a rising tail at the small size
% end of PSD. This originates from multiple scattering, and may not be
% real. In most cases, simply discard the smallest size bins. Note,
% particles below measurement range of instruments can produce this tail
% even when transmission is well above 30%. If transmission > 99%, your PSD
% will be noisy as there is not sufficient scattered energy on the
% detectors in the instrument.
% PARTICLE SIZE DISTRIBUTION: [LISST-100X only]  If your PSD shows a rising
% tail, and water depth is low, the scattering seen by ring detectors may
% be contaminated by ambient light, thus inventing fine particles. This
% problem does not exist in the LISST-200X; it is programmed to reject
% ambient light. Second, if the PSD has a spike at a single size, it is
% almost certainly an error; such spikes cannot arise from measured light
% scattering.
% NET SCATTERING ON RINGS: The most rigorous test of data quality, done
% only at the factory, is to test for smoothness of the net light on rings.
% Users who employ our Matlab software for data processing can view this
% property. Smooth net scattering can only arise from high quality data.
% The test is visual, so a little subjective. A saw-tooth pattern on net
% scattering can result from occasional fibrous material in the laser beam.
% In general, if data are saved as average of many samples, the sawtooth
% structure is minimal. Figure below shows many scans of high quality net
% particle scattering data.

% <http://www.sequoiasci.com/article/random-shaped-particles-lissts/>
% This new matrix was then compared to the matrix based on spherical
% particles.  The difference for large particles (> 16 µm) was minimal.
% However, for particles smaller than 16 µm, random-shaped particles
% scattered more light at larger angles than spheres of similar size.  This
% increased scattering at large angles became larger as particle size
% decreased.
% This phenomenon explains an often-observed rising tail in the fine end of
% the size distributions.  Fine, random-shaped particles scatter more light
% at large angles than spheres of similar size, so when the scattering
% pattern is inverted using a Kv matrix based on sphered the result is that
% fine particles are being ‘invented’.  This artificially causes an
% increase in particle volume in the fine end of the size distributions.
% Using the natural particle matrix to invert the scattering pattern causes
% the rising tail of fine particles to disappear.

% <http://www.sequoiasci.com/article/how-accurate-is-the-concentration-measurement-by-lisst-instruments/>

% <https://www.sequoiasci.com/article/the-influence-of-particles-outside-the-size-range-of-the-lisst/>
