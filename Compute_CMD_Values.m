%%
% user values
wavelength = 940;
cmd_value = .2;


%import the fit parameters from the latest calibration
close all;
folder = 'data';
files = dir(['./' folder '/*_fit_vals.csv']);

%find the most recent calibration

[~,idx] = sort([files.datenum], 'descend');
data = readmatrix(['./' folder '/' files(idx(1)).name]);

wavelengths = data(:,1);
photons_per_milliwatt = data(:,2);
a = data(:,3);
b = data(:,4);
c = data(:,5);
d = data(:,6);

f = @(a,b,c,d,x) a*sin(b*x+c)+d;
unique_wavelengths = unique(wavelengths);

this_wavelength_index = find(wavelengths == wavelength);

if(length(this_wavelength_index) > 1)
   error('Too many wavelength fit params that match in the database; possibly corrupted'); 
end

if(length(this_wavelength_index) < 1)
   error('There is no fit for that wavelength in the calibration'); 
end

photon_flux_at_wavelength = f(a(this_wavelength_index), b(this_wavelength_index), c(this_wavelength_index), d(this_wavelength_index), cmd_value);

%loop over the other values and solve for htis new peak value at lowest
%power wavelength
figure; hold on;
plot([0 1], [photon_flux_at_wavelength photon_flux_at_wavelength], 'r--');
y = photon_flux_at_wavelength;
for i = 1:length(unique_wavelengths)
    %plot the fit lines
    plot(0:0.0001:1, f(a(i),b(i),c(i),d(i),0:0.0001:1),'k');
    
    
    if photon_flux_at_wavelength > max(f(a(i),b(i),c(i),d(i),0:0.0001:1))
        equi_photon_cmd_value(i) = NaN;
    else
        x = (asin((y-d(i))/a(i))-c(i))/b(i);

        if(x < 0)
           x = x + (2*pi)/b(i); 
        end
        equi_photon_cmd_value(i) = x;
    end
    
    milliwatts_at_command_at_wavelength(i) = photon_flux_at_wavelength/photons_per_milliwatt(i);
    
    plot(equi_photon_cmd_value(i), photon_flux_at_wavelength, 'rx');
end
grid on;

figure; plot(unique_wavelengths,equi_photon_cmd_value, 'rx');
xlim([min(wavelengths)-10, max(wavelengths)+10]);
ylim([-0.1, 1.1]);
grid on;
figure; plot(unique_wavelengths,milliwatts_at_command_at_wavelength, 'bo');
grid on;