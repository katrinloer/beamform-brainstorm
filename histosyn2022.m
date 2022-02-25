%% Histogram from beamforming of synthetic wavefields
% Compute a histogram (velocity vs azimuth) from beamforming of synthetic,
% noise-like wavefields recorded on an array of choice and estimate the 
% apparent anisotropy by fitting a curve to the histogram
%--------------------------------------------------------------------------
% Katrin Loer
% Department of Geology and Geophysics
% University of Aberdeen
% United Kingdom
%
% Feb 2017
% last modified: Jan 2022
%--------------------------------------------------------------------------
% This script models different wave fields for a given number of 
% time windows, computes the corresponding beamresponse R and picks the
% maximum of R. All maxima for all time windows are then plotted in a
% histogram (velocity vs azimuth) and a curve is fitted to the distribution 
% (as in "real" beamforming analysis). Each histogram corresponds to a 
% given frequency.
% Wavefields are modelled in the frequency domain in terms of spectral
% amplitude and phase associated with the chosen wavefield parameters
% (velocity, frequency, azimuth) at each station of the array. Only the 
% vertical component is considered.
% Each wavefield can consist of multiple superimposed waves, which have the
% same wavenumber but different azimuths. Velocities are isotropic.  
% The resulting curve should give an idea about the apparent anisotropy
% intruduced as a result of the array geometry and the source distribution
% (that is, the superposition and azimuthal distribution of recorded
% waves)
%--------------------------------------------------------------------------
% Additional functions required:
% - extrema2.m
% - export_fig.m
%--------------------------------------------------------------------------
% Parameters to be set:
% Path to additional functions (l. 50)      - default is './'
% saveflag (l. 53)                          - default: true
% arrayflag (l. 67)                         - default: rand
% Array response parameters (ll. 164 ff.)   - default values see below
% Number of time windows nwin (l. 188)      - deafult: nwin=5000
% Wavefield design (ll. 192 ff.)            - default values see below
% Number of bootstrap resamples B (l. 322)  - default B=1000
% Plotting parameters (ll. 370 ff.)         - default values see below
% Path to figures/data (ll. 123+405)        - default is './'
%--------------------------------------------------------------------------

clear
close all

% addpath('/Users/s06kl9/Projects/FKanalysis/FK3C') % path to extrema2.m
% addpath('/Users/s06kl9/Projects/FKanalysis/SCRIPTS/altmany-export_fig-f0af704') % path to export_fig.m

saveflag = true;

tic % measures run time (for default values 6-7 min)

%% Load array coordinates
% Station coordinates must be in m (not lat/lon) and
% coords must be n x 2 double (n = number of stations)
% Choose/modify one of the predefined arrays:
% 'circ' = circular array
% 'rect' = rectangular array
% 'rand' = random array
% 'line' = linear array
% or load you own (here 'gemex')

arrayflag = 'rand';

switch arrayflag
    case 'gemex'
        % GEMEX
        load /Users/s06kl9/Projects/FKanalysis/GEMEX/coords_DB_centred.mat
    case 'circ'
        % Circular
        cth = (0:45:345)*pi/180;
        cr = 2000:2000:4000;
        X = kron(cr,cos(cth));
        Y = kron(cr,sin(cth));
        coords = [X(:) Y(:)];
    case 'rect'
        % Rectangular
        xmax = 7000;
        ymax = 1000;
        xstep = 2000;
        ystep = 2000;
        [X, Y] = meshgrid(-xmax:xstep:xmax,-ymax:ystep:ymax);
        coords = [X(:) Y(:)];
    case 'rand'
        % Random
        n = 16;
        smax = 12000;
        coords = smax*rand(n,2);
    case 'line'
        % Linear array
        xline = linspace(-2000,8000,8);
        yline = xline;
        ctest(:,1) = xline;
        ctest(:,2) = yline(3);
        ctest2(:,2) = yline;
        ctest2(:,1) = xline(3);
        ctest(3,:)=[];
        coords = cat(1,ctest,ctest2);
end
n = size(coords,1);
centre = mean(coords); % used below as reference point for source coordinates

% Plot array geometry
figure('Color','w'); 
plot(coords(:,1),coords(:,2),'bv','MarkerFaceColor','b','MarkerSize',12); 
axis equal
xmax = max(abs(coords(:,1))) + 500;
ymax = max(abs(coords(:,2))) + 500;
xlim([-xmax xmax])
ylim([-ymax ymax])
grid on
grid minor
axis equal
title('Array geometry')
set(gca,'Fontsize',18)
xlabel('x in m')
ylabel('y in m')

% Save array figure
if saveflag
    arrayfig = sprintf('./array_%s.png',arrayflag);
    export_fig(arrayfig,'-png','-r200')
end

%% Station pair orientation (used for plotting only)

kk = 0;
azi = 0;
dmax = 20000;
AA = zeros(length(coords),length(coords)); % matrix that stores station pair oreintations

for i = 1:length(coords)
    for j = i:length(coords)        
        if i ~= j
            x = coords(i,1) - coords(j,1);
            y = coords(i,2) - coords(j,2);
            % Distance between two stations
            dist = sqrt(x^2+y^2);
            if dist < dmax
                kk = kk+1;
                % Orientation of station pair
                azi_help = atan(y/x);
                if azi_help < 0
                    azi(kk) = pi + azi_help;
                else
                    azi(kk) = azi_help;
                end
                AA(i,j) = azi(kk);
                AA(j,i) = azi(kk)+pi;
            else
                continue
            end
        else
            AAmax = 9999;
            AA(i,j) = AAmax; % set diagonal (same station) to max to exclude from histogram later
        end
    end
end

%% Array response vector
% used to compute beamresponse of the array

kmax = 0.3 * 1/1000;        % maximum wavenumber
res = 201;                  % wavenumber resolution
kr = linspace(0,kmax,res)'; % range of wavenumbers
kth = (5:5:360)*pi/180;     % range of azimuths
k = kron(kr, [cos(kth(:)) sin(kth(:))] ); % complete range of wavenumber vectors

ak = 1/sqrt(n) * exp(1i*2*pi*coords*k'); % array response vector (requires wavenumber vectors and station coordinates)

%% "Synthetic" wavefield
% Define parameters in the "Wavefield design" box or use default
% - f and v: set according to observations you've made so far
% - d: for each time window, the dominant direction of the wavefield will
% be taken as a random number from the distribution defined in d; 
% by default this is the range [1 360] degree - you could limit that range,  
% for example, when you observed a dominant direction in your data
% - w: standard deviation of the distribution of these waves in degree
% - nk: number of waves superimposed in one time window
% NOTE: all waves will have the same wavenumber as defined by
% frequency and velocity, only their azimuths vary according to a normal
% distribution defined by d, w and nk

nwin = 5000;    % number of time windows: you could play with this number 
% to explore its influence and to find out if there is a "saturation point",
% i.e., a number beyond which the fitted curve doesn't change anymore

% Wavefield design
%--------------------------------------------------------------------------
% In the "time window loop" below, we model one wavefield for each time 
% window, computes the beamresponse and pick the maximum.
% One wavefield can consist of multiple superimposed waves from different
% directions; you can change the properties of the wavefield by varying the
% parameters in this "Wavefield design" box
% - fieldflag: chose 'rand' for random source locations or 'domi' for a 
% dominant direction of sources
% - d: dominant direction of sources
% - w: standard deviation of the distribution of these waves
% - nk: number of waves superimposed in one time window

fieldflag = 'rand';
f = 0.3;        % frequency in Hz
v = 2600;       % velocity in m/s
kr_s = f/v;     % horizontal wavenumber

switch fieldflag
    case 'domi'
        d = round(360 * rand(1,nwin)); % dominant direction of propagation (= mean of distribution)
        w = 45;         % standard deviation of distribution
        nkmax = 20;     % maximum number of events per time window
        nk = ceil(nkmax * rand(1,nwin)); % number of events per time window
    case 'rand'
        nkmax = 20;     % maximum number of events per time window
        nk = ceil(nkmax * rand(1,nwin)); % number of events per time window
end

ml = 5; % multiples of the wavelength, defines radius of source circle
rad_s = ml * 1/kr_s; % radius of cirlce
%--------------------------------------------------------------------------

%% Loop over time windows
% For each time window in the loop, we model the wavefield, compute the
% beamresponse and pick the maximum
% NOTHING TO BE CHANGED IN THE LOOP

kr_max = zeros(nwin,1);
kth_max = zeros(nwin,1);
for i = 1:nwin

    % 1) Model wavefield for ith time window
    %---------------------------------------
    switch fieldflag
        case 'domi'
            kth_s = ( d(i) + round( w * randn(1,nk)) ) * pi/180; % randn: normal distribution
            % kth_s = azimuthal distribution of superimposed waves using the mean
            % and standard deviation as defined above
        case 'rand'
            kth_s = round(360 * rand(1,nk(i))) * pi/180;
            % source distribution is random within 360 degrees
    end
    
    % Complete wavenumber vector
    k_s = kron(kr_s, [cos(kth_s(:)) sin(kth_s(:))]);
    
    % Source locations: The idea is to have the sources on a circle around 
    % the array, however, not exactly but slightly scattered to create a 
    % more realistic scenario
    coords_s = centre' + (rad_s + 1/kr_s*rand(1,nk(i))) .* [cos(pi+kth_s); sin(pi+kth_s)]; 
    
    % Wavefields in terms of spectral phases as a function of receiver
    % location relative to source location (coords-coords_s) and wave 
    % vector (k_s) (note that spectral amplitudes do not vary, i.e., they
    % are all equal to one so that all sources have the same contribution -
    % this could easily be modified though)
    s = zeros(n,nk(i));
    for kk = 1:nk(i)
        s(:,kk) = exp(1i*2*pi*(coords-coords_s(:,kk)')*k_s(kk,:)');
    end
    
    % Sum of all wavefields (coming from different azimuths)
    s = sum(s,2);
    
    % Normalize amplitudes for all stations
    s = s./abs(s);
    
    % 2) Compute beam response
    %-------------------------
    % Beam response after Riahi et al. (2013) eq. 3
    % R(k) = a(k)* (s s*) a(k), where (s s*) = S (spectral density matrix)
    R = ak.' * s; % gives direction of origin (ak'*s gives direction of propagation)
    R = R .* conj(R);
    
    R = reshape(R,size(kth,2),size(kr,1));
    
    % 3) Pick maximum
    %----------------
    [xmax,imax,~,~] = extrema2(R);
    [ii,jj] = ind2sub(size(R),imax(1)); % get first maximum
    
    kr_max(i) = kr(jj);
    kth_max(i) = kth(ii);
    
    if mod(i,100)==0
        fprintf('Done time window %d out of %d\n',i,nwin)
    end
    
end

% Convert from wave number to velocity
vr_max = f ./ kr_max;
vr = flipud( f ./ kr(kr > 0) );

% Convert from m/s to km/s
vr = vr ./ 1000;
vr_max = vr_max ./ 1000;

% Convert from radians to degree
kth = kth*180/pi;
kth_max = kth_max*180/pi;

%% Define anisotropy curve
% This curve will be fitted to the velocity distribution in the histogram

fit = @(b,x) b(1) + b(2)*cos(2*x*pi/180) + b(3)*sin(2*x*pi/180) + b(4)*cos(4*x*pi/180) +  b(5)*sin(4*x*pi/180);
% after Smith & Dahlen (1973)

%% Bootstrap
% Bootstrapping is any test or metric that uses random sampling with 
% replacement (e.g. mimicking the sampling process), and falls under the 
% broader class of resampling methods. Bootstrapping assigns measures of 
% accuracy (bias, variance, confidence intervals, prediction error, etc.) 
% to sample estimates.This technique allows estimation of the 
% sampling distribution of almost any statistic using random sampling 
% methods. (https://en.wikipedia.org/wiki/Bootstrapping_(statistics))

fprintf('Now bootstrapping...\n');

B = 1000; % number for bootstrap resamples

% Initial values and options for curve fitting:
X0 = zeros(1,5);
options = optimset('MaxFunEvals',400*length(X0),'MaxIter',400*length(X0));

i = 1;
while i <= B
    
    xyr = [kth_max, vr_max];
    N = length(vr_max); % number of samples used for resampling
    xyr_boot = datasample(xyr,N); % resampled dataset
    xr = xyr_boot(:,1);
    yr = xyr_boot(:,2);
    X0(1) = mean(yr); % initial value taken as average velocity
    
    % Define a function that computes the misfit between the measured
    % velocities (yr) and those computed from the anisotropy equation
    % (fit(b,xr)) using the L1-norm:
    fcn = @(b) sum(abs(fit(b,xr) - yr));
    
    % Find the anisotropy parameters (b) that minimise the misfit:
    [b,~,exitflag] = fminsearch(fcn,X0,options);
    
    % Store the anisotropy parameters for each bootstrap resample:
    if exitflag == 1
        asynth(i,:) = b;
        i = i + 1;
    end
end

% Mean anisotropy curve from all bootstrap resamples
y = fit(mean(asynth),kth);

% Amplitude of anisotropy
[amax, iamax] = max(y);
[amin, iamin] = min(y);
amean = mean(y);
amag = 1/2 * (amax-amin) / amean;
fprintf('Magnitude of apparent anisotropy: %1.1f per cent\n',amag*100)

% Fast direction of anisotropy
afast = kth(iamax);
if afast > 180
    afast = afast - 180;
end
fprintf('Fast direction of apparent anisotropy: %3.0f degree\n',afast)

%% Plot

% Histogram
figure('Color','white');
hold on;
set(gcf,'Position',[100 100 1000 500])
histogram2(kth_max,vr_max,kth,vr,'DisplayStyle','tile','ShowEmptyBins','On');
ylim([min(vr) 4]) % adapt y-axis as you see fit for your velocity range
xlabel('azimuth in degree')
ylabel('velocity in km/s')
color_love = colorbar;
color_love.Label.String = 'number of detections';
set(color_love,'Fontsize',32)
set(gca,'Fontsize',24,'XTick',45:45:360)
set(gca,'YTick',1:0.5:4)

% Plot isotropic (expected) velocity as black line
line([5 360],[v/1000 v/1000],'Color','k')

% Plot curves fitted to bootstrap resamples as grey curves
for i = 1:B
    plot(kth,fit(asynth(i,:),kth),'LineWidth',1,'color',[0.5,0.5,0.5]);
end

% Plot mean of bootstrap curves
plot(kth,fit(mean(asynth),kth),'w','LineWidth',2);

% Plot station pair orientation
yyaxis right
set(gca,'YColor','k')
histogram(AA*180/pi,'BinLimits',[0 360],'BinWidth',10,'FaceColor',[0.5,0.5,0.5])
ylim([0 n*4])
ylabel('number of station pairs')
xlim([min(kth) max(kth)])

%% Save data & figure

if saveflag

    filename = sprintf('./asynth_%s_f%1.2f_nwin%d_nkrand.mat',arrayflag,f,nwin);
    save(filename,'asynth','f','v','nk')
    
    figname = sprintf('./histogram_synth_%s_f%1.2f_v%04d_nkrand_nwin%05d.png',arrayflag,f,v,nwin);
    export_fig(figname,'-png','-r200')

end

%% EOF
toc
