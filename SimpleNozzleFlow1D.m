function [Flow_parameters] = SimpleNozzleFlow1D(Ambient_Input, Nozzle_Input, Fuel_Input, chamber_pressure, Booleans, n, plot_points_number)

%Simulation_Results.pressure(n) - chamber pressure
%Fuel_Input.temperature_combustion - chamber temperature
chamber_density        = chamber_pressure / (Fuel_Input.gas_constant * Fuel_Input.temperature_combustion);
chamber_sonic_velocity = sqrt(Fuel_Input.temperature_combustion * Fuel_Input.kappa * Fuel_Input.gas_constant);

Flow_parameters.exit_mach = MachBisection(Fuel_Input.kappa, Nozzle_Input.exit_diameter^2 / Nozzle_Input.throat_diameter^2, 0.00001, [1 15]);
Flow_parameters.exit_pressure = chamber_pressure * PressureFromMachRatio(Flow_parameters.exit_mach, Fuel_Input.kappa);

if true   %"Normal" flow
    Flow_parameters.x_convergent = linspace(0, Nozzle_Input.convergent_length, plot_points_number);%Nozzle shape functions
    Flow_parameters.x_divergent  = linspace(Nozzle_Input.convergent_length + (Nozzle_Input.convergent_length + Nozzle_Input.divergent_length) / plot_points_number, Nozzle_Input.convergent_length + Nozzle_Input.divergent_length, plot_points_number);
    
    Flow_parameters.d_convergent = interp1(Nozzle_Input.points(:,1), Nozzle_Input.points(:,2) .* 2, Flow_parameters.x_convergent, 'linear', 'extrap');
    Flow_parameters.d_divergent  = interp1(Nozzle_Input.points(:,1), Nozzle_Input.points(:,2) .* 2, Flow_parameters.x_divergent, 'linear', 'extrap');

    Flow_parameters.entry_mach              = MachBisection(Fuel_Input.kappa, Nozzle_Input.entry_diameter^2 / Nozzle_Input.throat_diameter^2, 0.01, [eps 1]);
    Flow_parameters.mach_numbers_convergent = linspace(Flow_parameters.entry_mach,1,plot_points_number);
    
    convergent_area_ratios_from_mach        = AreaRatioFromMach(Fuel_Input.kappa ,Flow_parameters.mach_numbers_convergent);
    Flow_parameters.mach_numbers_convergent = interp1(convergent_area_ratios_from_mach,Flow_parameters.mach_numbers_convergent,(Flow_parameters.d_convergent ./ Nozzle_Input.throat_diameter) .^2, "spline");
    Flow_parameters.mach_numbers_divergent  = linspace(1+(Flow_parameters.exit_mach-1)/plot_points_number,Flow_parameters.exit_mach,plot_points_number);
    
    divergent_area_ratios_from_mach         = AreaRatioFromMach(Fuel_Input.kappa ,Flow_parameters.mach_numbers_divergent);
    Flow_parameters.mach_numbers_divergent  = interp1(divergent_area_ratios_from_mach, Flow_parameters.mach_numbers_divergent,(Flow_parameters.d_divergent ./ Nozzle_Input.throat_diameter) .^2, "spline");
    
    %pressures
    Flow_parameters.pressures_convergent = chamber_pressure .* PressureFromMachRatio(Flow_parameters.mach_numbers_convergent, Fuel_Input.kappa);
    Flow_parameters.pressures_divergent  = chamber_pressure .* PressureFromMachRatio(Flow_parameters.mach_numbers_divergent, Fuel_Input.kappa);
    Flow_parameters.throat_pressure      = Flow_parameters.pressures_convergent(plot_points_number);
    
    %temperatures
    Flow_parameters.temperatures_convergent = Fuel_Input.temperature_combustion .* TemperatureFromMachRatio(Flow_parameters.mach_numbers_convergent, Fuel_Input.kappa);
    Flow_parameters.temperatures_divergent  = Fuel_Input.temperature_combustion .* TemperatureFromMachRatio(Flow_parameters.mach_numbers_divergent, Fuel_Input.kappa);
    Flow_parameters.throat_temperature      = Flow_parameters.temperatures_convergent(plot_points_number);
    Flow_parameters.exit_temperature        = Flow_parameters.temperatures_divergent(plot_points_number);
    
    %densities
    Flow_parameters.densities_convergent = chamber_density .* DensityFromMachRatio(Flow_parameters.mach_numbers_convergent, Fuel_Input.kappa);
    Flow_parameters.densities_divergent  = chamber_density .* DensityFromMachRatio(Flow_parameters.mach_numbers_divergent, Fuel_Input.kappa);
    
    %velocities
    Flow_parameters.velocities_convergent = Flow_parameters.mach_numbers_convergent .* chamber_sonic_velocity .* SonicVelocityFromMachRatio(Flow_parameters.mach_numbers_convergent, Fuel_Input.kappa);
    Flow_parameters.velocities_divergent  = Flow_parameters.mach_numbers_divergent .* chamber_sonic_velocity .* SonicVelocityFromMachRatio(Flow_parameters.mach_numbers_divergent, Fuel_Input.kappa);
    Flow_parameters.throat_velocity       = Flow_parameters.velocities_convergent(plot_points_number);
    Flow_parameters.exit_velocity         = Flow_parameters.velocities_divergent(plot_points_number);
    
    if (~isreal(Flow_parameters.exit_velocity) && Booleans.running_nozzle_optimization == false)
        warning('Complex nozzle_velocity in SimpleNozzleFlow1D.m. n == %d',n);
        Flow_parameters.exit_velocity = 0;
    end 
    
    [Flow_parameters.thrust, Flow_parameters.mass_flow] = exit(Flow_parameters.densities_convergent(plot_points_number), Flow_parameters.throat_velocity, Flow_parameters.exit_velocity, Nozzle_Input, Flow_parameters.exit_pressure, Ambient_Input, Booleans);
    Flow_parameters.throat_mach = 1;
end

end

%% Functions

% exit
function [thrust, mass_flow] = exit(throat_density, throat_velocity, exit_velocity, Nozzle_Input, exit_pressure, Ambient_Input, Booleans)
    mass_flow = throat_density * throat_velocity * Nozzle_Input.throat_diameter^2 * pi / 4;
    thrust    = exit_velocity * mass_flow;
   
    if Booleans.correction_factors
        dl = Nozzle_Input.points(end,1) - Nozzle_Input.points(end - 1,1);
        dd = Nozzle_Input.points(end,2) - Nozzle_Input.points(end - 1,2);
        thrust = thrust * (1 + cos(atan(dd / dl))) / 2;
    end

    if Booleans.thrust_pressure_term
        thrust = thrust + (exit_pressure - Ambient_Input.pressure) * (Nozzle_Input.exit_diameter ^ 2) * pi / 4;
        if thrust < 0
            thrust = eps;
        end
    end
    thrust    = thrust * Nozzle_Input.efficiency;
end

% MachBisection
function [value] = MachBisection(kappa, area_ratio, errorr, values)
    fun = @(M) (2 ./ (kappa + 1)).^((kappa + 1) ./ (2 * (kappa - 1))) .* (1 + ((kappa - 1) ./ 2 * M.^2)).^((kappa + 1) ./ (2 * (kappa - 1))) ./ M - area_ratio;
    fv(1) = fun(values(1));
    fv(2) = fun(values(2));
    
    if fv(1) ~= 0 
        if fv(2) ~= 0
            values(3) = (values(1) + values(2)) / 2;
            fv(3)     = fun((values(1) + values(2)) / 2);
            
            while ((abs(fv(3)) > errorr))
                values(3) = (values(1) + values(2)) / 2;
                fv(3)     = fun((values(1) + values(2)) / 2);
                x         = values(3);

                if fv(3) * fv(2) < 0
                    values(1) = x;
                    fv(1) = fun(values(1));
                    fv(2) = fun(values(2));
                    fv(3) = fun((values(1) + values(2)) / 2);

                elseif fv(3) * fv(1) < 0
                    values(2) = x;
                    fv(1) = fun(values(1));
                    fv(2) = fun(values(2));
                    fv(3) = fun((values(1) + values(2)) / 2);

                else
                    value = 1;
                    if abs(fun(value)) < 100 * eps
                        % warning("Bisection failed but fun(1) is close to 0 so we assume that 1 is ok");
                        return;

                    else
                        error("Bisection_failed");
                    end
                end
            end
            
            value = values(3);

        else
            value = values(2);
        end

    else
        value = values(1);
    end

end

%% Isentropic relations

% AreaRatioFromMach
function [area_ratio] = AreaRatioFromMach(kappa, mach)
    area_ratio = (2 ./ (kappa + 1)).^((kappa + 1) / (2 * (kappa - 1))) .* (1 + ((kappa - 1) ./ 2 .* mach.^2)).^...
     ((kappa + 1) ./ (2 .* (kappa - 1))) ./ mach;
end

% TemperatureFromMachRatio
function [ratio] = TemperatureFromMachRatio(M, kappa)
    ratio = 1 ./ (1 + ((kappa - 1) / 2) .* M.^2);
end

% PressureFromMachRatio
function [ratio] = PressureFromMachRatio(M, kappa)
    ratio = (1 + ((kappa - 1) / 2) .* M.^2) .^ (kappa / (1 - kappa));
end

% DensityFromMachRatio
function [ratio] = DensityFromMachRatio(M, kappa)
    ratio = (1 + ((kappa - 1) / 2) .* M.^2) .^ (1 / (1 - kappa));
end

% SonicVelocityFromMachRatio
function [ratio] = SonicVelocityFromMachRatio(M, kappa)
    ratio = 1 ./ sqrt(1 + ((kappa - 1) / 2) * M.^2);
end