%% FUNCTION find_turning_points
  function TP = find_turning_points(time,depth)
    dt = 5./24./60./60; % 5-second bin size
    % bin the pressure
    tbin = [nanmin(time):dt:nanmax(time)]';
    dbin = binaverage( time, depth, tbin );
    % interpolate nans
    dbin = fillInvalidValues(dbin,'linear');
    
    % use the binned data to calculate first and second derivatives of pressure
    dpdt = gradient( dbin, tbin );   % vertical velocity
    d2pdt2 = gradient( dpdt, tbin ); % double derivative of pressure
    
    % make an up-down vector from the binned pressure data
    % (+1 rising, -1 sinking, 0 constant depth)
    udvec = nan*dbin;
    udvec( sign( dpdt ) <0 ) = +1; % find all the ascending data points
    udvec( sign( dpdt ) >0 ) = -1; % find all the descending data points
    udvec(isnan(udvec)) = 0;
    % curvature;
    k = abs(d2pdt2) ./ (1 + dpdt.^2).^1.5;
    
    % find all the turning points in the up-down vector
    tpb = find( diff( udvec) ~=0 );
    %tpb = find( abs(diff( udvec)) == 2);
    % plot( tbin(tpb), dbin(tpb), 'k*','MarkerSize',10); hold on
    % shift one indice to the right
    tpb = tpb+1;
    % now to only keep the turning points that define the start of a new
    % profile. Empirically decided that profiles with less than 5 data points
    % are small 'jogs' undergone by the glider rather than a total profile.
    faketurns = (diff( tpb) < 5); % finplotd fake turns dur to jogs
    tpb(faketurns) = nan;  % throw out fake turns.
    tpb = tpb(~isnan( tpb ));
    
    % Define whether each profile is up or down
    updowns = udvec( tpb );
    % flag bad indices where profiles are same direction
    badinds = find( diff(updowns) ==0 );    
    badinds = badinds+1;
    
    % now remove the bad indices
    updowns(badinds) = nan;
    updowns = updowns( ~isnan( updowns ));
    tpb(badinds) = nan;
    tpb = tpb( ~isnan( tpb ));
    % keyboard
    %  now use find the turning points in the original data vector
    TP = [];
    TP = interp1( time, 1:length( time ), tbin(tpb) );
    TP = floor( TP );
  end%% FUNCTION find_turning_points
