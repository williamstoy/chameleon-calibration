function [photon_flux_at_wavelength, equi_photon_cmd_value, unique_wavelengths, milliwatts_at_command_at_wavelength, bias_voltages, front_panel_bias, filename] = Compute_CMD_Values_Func(wavelength, cmd_value, hSinFitPlot, hEquiphotonFluxPlot, hPowerAtEquiphotonFlux)
    cla(hSinFitPlot);
    cla(hEquiphotonFluxPlot);
    cla(hPowerAtEquiphotonFlux);

    %import the fit parameters from the latest calibration
    folder = 'data';
    fit_values_files = dir(['./' folder '/*_fit_vals.csv']);
    bias_min_files = dir(['./' folder '/chameleon_bias_minimum_sweep_raw_data*.csv']);

    %find the most recent fit calibration
    [~,idx] = sort([fit_values_files.datenum], 'descend');
    filename = fit_values_files(idx(1)).name;
    fit_data = readmatrix(['./' folder '/' filename]);
    
    %find the most recent bias calibration
    [~,idx] = sort([bias_min_files.datenum], 'descend');
    filename = bias_min_files(idx(1)).name;
    bias_data = readmatrix(['./' folder '/' filename]);
    
    bias_wavelength = bias_data(:,1);
    front_panel_bias = bias_data(:,3);
    bias_voltages = bias_data(:,4);
    
    [~, bias_indicies] = sort(bias_wavelength);
    front_panel_bias = front_panel_bias(bias_indicies);
    bias_voltages = bias_voltages(bias_indicies);

    % process the fit data
    wavelengths = fit_data(:,1);
    photons_per_milliwatt = fit_data(:,2);
    a = fit_data(:,3);
    b = fit_data(:,4);
    c = fit_data(:,5);
    d = fit_data(:,6);

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
    
    plot(hSinFitPlot, [0 1], [photon_flux_at_wavelength photon_flux_at_wavelength], 'r--'); hold(hSinFitPlot, 'on'); 
    y = photon_flux_at_wavelength;
    for i = 1:length(unique_wavelengths)
        %plot the fit lines
        if(unique_wavelengths(i) == wavelength)
            plot_color = 'r';
            line_width = 3;
        else
            plot_color = 'k';
            line_width = 1;
        end
        plot(hSinFitPlot, 0:0.0001:1, f(a(i),b(i),c(i),d(i),0:0.0001:1),plot_color, 'LineWidth', line_width);

        if photon_flux_at_wavelength > max(f(a(i),b(i),c(i),d(i),0:0.0001:1))
            equi_photon_cmd_value(i) = NaN;
            milliwatts_at_command_at_wavelength(i) = NaN;
        else
            x = (asin((y-d(i))/a(i))-c(i))/b(i);

            if(x < 0)
               x = x + (2*pi)/b(i); 
            end
            equi_photon_cmd_value(i) = x;
            milliwatts_at_command_at_wavelength(i) = photon_flux_at_wavelength/photons_per_milliwatt(i);
        end

        plot(hSinFitPlot, equi_photon_cmd_value(i), photon_flux_at_wavelength, 'rx', 'MarkerSize', 10);
    end
    plot(hSinFitPlot, [cmd_value, cmd_value], [0, photon_flux_at_wavelength], 'r--');
    xlim(hSinFitPlot, [0,1]);
    ylim(hSinFitPlot, [0,inf]);
    grid on;

    
    plot(hEquiphotonFluxPlot, unique_wavelengths,equi_photon_cmd_value, 'rx', 'MarkerSize', 10); hold(hEquiphotonFluxPlot, 'on');
    equi_photon_cmd_value_at_input_wavelength = equi_photon_cmd_value(unique_wavelengths == wavelength);
    plot(hEquiphotonFluxPlot, [wavelength, wavelength], [-0.1, equi_photon_cmd_value_at_input_wavelength], 'r--');
    plot(hEquiphotonFluxPlot, [min(wavelengths)-10, wavelength], [equi_photon_cmd_value_at_input_wavelength, equi_photon_cmd_value_at_input_wavelength], 'r--');
    xlim(hEquiphotonFluxPlot, [min(wavelengths)-10, max(wavelengths)+10]);
    ylim(hEquiphotonFluxPlot, [-0.1, 1.1]);
    grid(hEquiphotonFluxPlot, 'on');
    
    plot(hPowerAtEquiphotonFlux, unique_wavelengths,milliwatts_at_command_at_wavelength, 'bo', 'MarkerSize', 10); hold(hPowerAtEquiphotonFlux, 'on'); 
    milliwatts_at_command_at_input_wavelength = milliwatts_at_command_at_wavelength(unique_wavelengths == wavelength);
    plot(hPowerAtEquiphotonFlux, [wavelength, wavelength], [0, milliwatts_at_command_at_input_wavelength], 'r--');
    plot(hPowerAtEquiphotonFlux, [min(wavelengths)-10, wavelength], [milliwatts_at_command_at_input_wavelength, milliwatts_at_command_at_input_wavelength], 'r--');
    grid(hPowerAtEquiphotonFlux, 'on');
    xlim(hPowerAtEquiphotonFlux, [min(wavelengths)-10, max(wavelengths)+10]);
    ylim(hPowerAtEquiphotonFlux, [0, 1000]);
    
    milliwatts_at_command_at_wavelength = round(abs(milliwatts_at_command_at_wavelength'));
    photon_flux_at_wavelength = abs(photon_flux_at_wavelength);
    equi_photon_cmd_value = round(abs(equi_photon_cmd_value'),3);
end