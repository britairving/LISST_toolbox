function ret = LISST_bin_sizes(type2)
%FUNCTION LISST_bin_sizes
%
%  Syntax:
%    bins = LISST_bin_sizes(opt.inst.type2)
%
%     where type2 is
%       1 for type A (discontinued)        | type A (5-500 µm)
%       2 for type B                       | type B (1.25-250 µm size range) 
%       3 for type C                       | type C (2.5-500 µm size range) 
%       4 for FLOC (discontinued)          | FLOC   (7.5-1500 µm size range)
%       21 for randomly shaped type B
%       31 for randomly shaped type C.
%
%  Description:
%    Modified from compute_mean.m Sequoia script
%
%  Refereces:
%    compute_mean.m
%
%  Authors:
%    Brita K Irving  <bkirving@alaska.edu>
%%
rho=200^(1/32);
switch type2 
  case 1 %% type A (5-500 µm)
    rho=100^(1/32);
    bins(:,1) = 5*rho.^([0:31]); %lower limit for type A
    bins(:,2) = 5*rho.^([1:32]); % upper limit for type A
    bins(:,3) = sqrt(bins(:,1).*bins(:,2));%mid-point for type A
    % dias32 = bins(:,3);%The midpoint of the size bins is being used for computation of means and stds
    % upperBins = bins(:,2);%Upper bin limits are being used for computing D50 (median) and other percentiles
    bin64 = 17;% this is the bin number containing 64µm particles. It is being used for computing the silt fraction later on.
    % return info
    ret.low_limit = bins(:,1);
    ret.upp_limit = bins(:,2);
    ret.mid_point = bins(:,3);
    ret.bin_size  = ret.upp_limit(1:end) - ret.low_limit(1:end);
    ret.bin64     = bin64;
  case 2 %% type B (1.25-250 µm size range) 
    bins(:,1) = 1.25*rho.^([0:31]); %lower limit for type B
    bins(:,2) = 1.25*rho.^([1:32]); % upper limit for type B
    bins(:,3) = sqrt(bins(:,1).*bins(:,2));%mid-point for type B
    % return info
    ret.low_limit = bins(:,1);
    ret.upp_limit = bins(:,2);
    ret.mid_point = bins(:,3);
    ret.bin_size  = ret.upp_limit(1:end) - ret.low_limit(1:end);
    ret.bin64     = 23;
  case 3 %% type C (2.5-500 µm size range) 
    bins(:,1) = 2.5*rho.^([0:31]); %lower limit for type C
    bins(:,2) = 2.5*rho.^([1:32]); % upper limit for type C
    bins(:,3) = sqrt(bins(:,1).*bins(:,2));%mid-point for type C
    % return info
    ret.low_limit = bins(:,1);
    ret.upp_limit = bins(:,2);
    ret.mid_point = bins(:,3);
    ret.bin_size  = ret.upp_limit(1:end) - ret.low_limit(1:end);
    ret.bin64     = 19;
  case 4 %% FLOC   (7.5-1500 µm size range)
    bins(:,1) = 7.5*rho.^([0:31]); %lower limit for type FLOC
    bins(:,2) = 7.5*rho.^([1:32]); %upper limit for type FLOC
    bins(:,3) = sqrt(bins(:,1).*bins(:,2));%mid-point for type FLOC
    % return info
    ret.low_limit = bins(:,1);
    ret.upp_limit = bins(:,2);
    ret.mid_point = bins(:,3);
    ret.bin_size  = ret.upp_limit(1:end) - ret.low_limit(1:end);
    ret.bin64     = 12;
  case 21 %% TYPE B randomnly shaped particles
    ret.ds=1*1.18.^(0:1:32); % from InLineAnalysis lib/processLISST.m and instrments/LISST.m
    dias32 = [1.0863095E+00  1.1800684E+00; 1.2819196E+00  1.3925615E+00; 1.5127528E+00  1.6433178E+00;...
      1.7851518E+00  1.9392274E+00; 2.1066013E+00  2.2884211E+00; 2.4859336E+00  2.7004934E+00;...
      2.9335718E+00  3.1867670E+00; 3.4618154E+00  3.7606031E+00; 4.0851790E+00  4.4377689E+00;...
      4.8207907E+00  5.2368710E+00; 5.6888629E+00  6.1798660E+00; 6.7132474E+00  7.2926647E+00;...
      7.9220913E+00  8.6058433E+00; 9.3486097E+00  1.0155484E+01; 1.1031999E+01  1.1984166E+01;...
      1.3018514E+01  1.4142136E+01; 1.5362737E+01  1.6688688E+01; 1.8129081E+01  1.9693793E+01;...
      2.1393555E+01  2.3240023E+01; 2.5245859E+01  2.7424818E+01; 2.9791841E+01  3.2363161E+01;...
      3.5156411E+01  3.8190744E+01; 4.1486970E+01  4.5067691E+01; 4.8957463E+01  5.3182959E+01;...
      5.7773156E+01  6.2759530E+01; 6.8176276E+01  7.4060540E+01; 8.0452671E+01  8.7396504E+01;...
      9.4939656E+01  1.0313385E+02; 1.1203529E+02  1.2170500E+02; 1.3220931E+02  1.4362023E+02;
      1.5601603E+02  1.6948170E+02; 1.8410959E+02  2.0000000E+02];%mid points (column 1) and upper bins (column 2) for type B, randomly shaped
    % return info
    ret.bin_size  = (ret.ds(2:end) - ret.ds(1:end-1)); % .* ones(size(dias32,1),1);
    ret.low_limit = [];
    ret.upp_limit = dias32(:,2);
    ret.mid_point = dias32(:,1);
    ret.bin64     = 25;
  case 31 %% TYPE C randomnly shaped particles
    ret.ds = [1.90, 2.25,2.65,3.13,3.69,4.35,5.14,6.06,7.15,8.44,9.96,11.8,13.9,...
          16.4,19.3,22.8,26.9,31.8,37.5,44.2,52.2,61.6,72.7,85.7,101,119,141,166,196,232,273,322,381];
    %mid points (column 1) and upper bins (column 2) for type C, randomly shaped
    dias32 = [2.05970200000000,2.24514560000000;2.43230210000000,2.64942550000000;2.87230560000000,3.12650330000000;...
      3.39190560000000,3.68948780000000;4.00550140000000,4.35384800000000;4.73009660000000,5.13783860000000;...
      5.58577110000000,6.06300100000000;6.59623700000000,7.15475600000000;7.78949630000000,8.44310160000000;...
      9.19861620000000,9.96343760000000;10.8626460000000,11.7575380000000;12.8276990000000,13.8746990000000;...
      15.1482290000000,16.3730940000000;17.8885440000000,19.3213720000000;21.1245810000000,22.8005400000000;...
      24.9460180000000,26.9061980000000;29.4587530000000,31.7511540000000;34.7878410000000,37.4685340000000;...
      41.0809620000000,44.2154340000000;48.5125080000000,52.1772370000000;57.2884200000000,61.5727100000000;...
      67.6518960000000,72.6600100000000;79.8901240000000,85.7437830000000;94.3422470000000,101.183530000000;...
      111.408760000000,119.403490000000;131.562600000000,140.904290000000;155.362280000000,166.276700000000;...
      183.467320000000,196.217880000000;216.656550000000,231.550520000000;255.849720000000,273.245460000000;...
      302.132940000000,322.448340000000;356.788790000000,380.511100000000];
    % return info
    ret.bin_size  = (ret.ds(2:end) - ret.ds(1:end-1)); % .* ones(size(dias32,1),1);
    ret.low_limit = [];
    ret.upp_limit = dias32(:,2);
    ret.mid_point = dias32(:,1);
    ret.bin64     = 21;
  otherwise
    error(['You must specify a number: 1-4, 21 or 31, NOT ',num2str(type),'.']);
end

end
