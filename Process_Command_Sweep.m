clc; clear; close all;
folder = 'data';
filename_no_filetype = 'chameleon_cmd_sweep_raw_data_2020-02-03_21-38';
data = readmatrix(['./' folder '/' filename_no_filetype '.csv'], 'NumHeaderLines',1);

wavelengths = data(:,1);
command = data(:,5);
photons_per_pulse = data(:,7);
power = data(:,6);

unique_wavelengths = unique(wavelengths);


f = @(a,b,c,d,x) a*sin(b*x+c)+d;

figure; hold on;
fit_params = [];
photons_per_milliwatt = [];
%% loop over the unique wavelengths and fit a sinusoid
for i = 1:length(unique_wavelengths)
   %clf;
   this_wavelength = unique_wavelengths(i);
   this_command = command(wavelengths == this_wavelength);
   this_photons_per_pulse = photons_per_pulse(wavelengths == this_wavelength);
   photons_per_milliwatt(i) = mean(this_photons_per_pulse./power(wavelengths == this_wavelength));
   
   amplitude = (max(this_photons_per_pulse)-min(this_photons_per_pulse));

   [~, max_location] = max(this_photons_per_pulse);
   [~, min_location] = min(this_photons_per_pulse);
   min_cmd = this_command(min_location);
   max_cmd = this_command(max_location);
   voltage_period = max_cmd - min_cmd;
   
   voltage_period = 1; %overwrite because we may not use the full range

   
   start_point = [amplitude/2*(1/max_cmd), 1.5*pi*voltage_period, 1.5*pi*voltage_period, amplitude/2*(1/max_cmd)];
   lower_bound = [0, 0.8*pi, 0.8*pi, 0];
   upper_bound = [inf, 1.8*pi, 1.8*pi, inf];
   this_fit = fit(this_command, this_photons_per_pulse, f, 'StartPoint', start_point, 'Lower', lower_bound, 'Upper', upper_bound);
   
   fit_params(i,:) = [this_fit.a, this_fit.b, this_fit.c, this_fit.d];
   
   plot(this_command, this_photons_per_pulse, 'ko'); hold on;
   plot(this_fit, 'k');
   %title([num2str(unique_wavelengths(i)) ' nm']);
   %pause;
end

writematrix([unique_wavelengths,photons_per_milliwatt',fit_params],['./' folder '/' filename_no_filetype '_fit_vals.csv']);
