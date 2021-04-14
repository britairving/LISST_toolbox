function [profile_index, profile_direction] = select_downcast(time,depth)
makefig;
profile_index     = nan(size(time));
profile_direction = nan(size(time));
done = 0;
while ~done
  plot( time, depth, 'k-'); hold on; axis ij; datetick
  title('Manually select profile limits')
  fprintf('Select profile starting place on plot\n')
  plot_select = ginput(1);
  startind = find(time >= plot_select(1),1);
  plot( time(startind), depth(startind), 'kp', 'markerfacecolor', 'g', 'markersize', 15)
  fprintf('Select profile ending place on plot\n')
  plot_select = ginput(1);
  endind = find(time <= plot_select(1),1,'last');
  plot( time(endind), depth(endind), 'kp', 'markerfacecolor', 'r', 'markersize', 15)
  profile_index(startind:endind) = 1;
  profile_direction(startind:endind) = -1;
  plot( time(startind:endind), depth(startind:endind), 'b-', 'LineWidth', 2)
  chc = input('  Is this correct? <1/0> ');
  if chc == 0
    done = 0;
    clf
  else
    close
    done = 1;
  end
end
end %% MAIN FUNCTION