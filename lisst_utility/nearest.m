function idx=nearest(dat,val)
 for ii = 1:length(val)
    idx(ii) = min(find(abs(dat-val(ii))==min(abs(dat-val(ii)))));
 end