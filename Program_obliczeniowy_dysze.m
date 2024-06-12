%% Inicjalizacja
[status, pythonPath] = system('where python');
if status == 0
    pythonPath = strtrim(pythonPath);
    pyenv('Version', pythonPath);
else
    error('Nie znaleziono pythona');
end
currentFolder = fileparts(mfilename('fullpath'));
if count(py.sys.path, currentFolder) == 0
    insert(py.sys.path, int32(0), currentFolder);
end
combustion_module = py.importlib.import_module('hydro_combustion');

%% Parametry wejściowe
pressure = 10000000;  %Ciśnienie komory [Pa]
hydrogen_fraction = 2/3;
oxygen_fraction = 1/3;

Nozzle_Input.entry_diameter = 0.1;  %średnica wejściowa dyszy [m]
Nozzle_Input.throat_diameter = 0.03;  %średnica krytyczna dyszy [m]
Nozzle_Input.throat_rounding_radius_1 = 0.01;    %Nie zmieniać
Nozzle_Input.throat_rounding_radius_2 = 0.02;    %Nie zmieniać
Nozzle_Input.exit_diameter = 0.06;  %średnica wyjściowa dyszy [m]
Nozzle_Input.convergent_length = 0.1;  %długość części zbieżnej dyszy [m]
Nozzle_Input.throat_length = 0.01;  %Długość odcinka o stałym przekroju w gardzieli [m]
Nozzle_Input.divergent_length = 0.15;  %śdługość części rozbieżnej dyszy [m]
Nozzle_Input.efficiency = 1;        %Wydajność dyszy

Booleans.correction_factors = false;    %Nie zmieniać
Booleans.thrust_pressure_term = true;    %Nie zmieniać
Booleans.enable_nozzle_erosion = false;    %Nie zmieniać

Ambient_Input.pressure = 100000;    %Ciśnienie otoczenia [Pa]
Ambient_Input.temperature = 300;    %Nie zmieniać

Nozzle_Input.shape_case = "conical";    %Nie zmieniać
Nozzle_Input = InitialiseNozzleShape(Nozzle_Input);

%% Wywołanie reakcji spalania
result = combustion_module.calculate_combustion_products(pressure, hydrogen_fraction, oxygen_fraction);
temperature = double(result{1});
pressure = double(result{2});
composition = result{3};
composition_keys = cell(py.list(composition.keys()));
composition_values = cell(py.list(composition.values()));
fprintf('Temperatura spalin: %.2f K\n', temperature);
fprintf('Ciśnienie spalin: %.2f Pa\n', pressure);
disp('Skład spalin:');
for i = 1:length(composition_keys)
    species = string(composition_keys{i});
    fraction = double(composition_values{i});
    fprintf('  %s: %.6f\n', species, fraction);
end
composition_frac = [];
for i = 1:length(composition_keys)
    composition_frac = [composition_frac, double(composition_values{i})];
end
gas_results_py = combustion_module.hydro_product_properties(temperature, pressure, composition_frac(1), composition_frac(2), composition_frac(3), composition_frac(4), composition_frac(5), composition_frac(6), composition_frac(7), composition_frac(8));
molecular_weight_py = gas_results_py{1};
gas_constant_py = gas_results_py{2};
gamma_py = gas_results_py{3};
viscosity_py = gas_results_py{4};
thermal_conductivity_py = gas_results_py{5};
molecular_weight = double(molecular_weight_py);
gas_constant = double(gas_constant_py);
gamma = double(gamma_py);
viscosity = double(viscosity_py);
thermal_conductivity = double(thermal_conductivity_py);
Fuel_Input.temperature_combustion = temperature;
Fuel_Input.kappa = gamma;
Fuel_Input.gas_constant = gas_constant;
%% Obliczenia przepływowe
figure("Position",[200,100,1100,600]);
tiledlayout(3,3);
nexttile;
plot(Nozzle_Input.points(:,1), - Nozzle_Input.points(:,2), "Color", [0 0 0]);
hold on;
plot(Nozzle_Input.points(:,1), Nozzle_Input.points(:,2), "Color", [0 0 0]);
grid minor;
axis equal;
hold off;
xlabel("distance from entry [m]");
ylabel("distance from main axis [m]");
title("Nozzle shape");

[Flow_parameters] = SimpleNozzleFlow1D(Ambient_Input, Nozzle_Input, Fuel_Input, pressure, Booleans, 1, 400);

thrust = Flow_parameters.thrust;

x_coordinates = [Flow_parameters.x_convergent, Flow_parameters.x_divergent];
diameters = [Flow_parameters.d_convergent, Flow_parameters.d_divergent];
pressures = [Flow_parameters.pressures_convergent, Flow_parameters.pressures_divergent];
temperatures = [Flow_parameters.temperatures_convergent, Flow_parameters.temperatures_divergent];
densities = [Flow_parameters.densities_convergent, Flow_parameters.densities_divergent];
velocities = [Flow_parameters.velocities_convergent, Flow_parameters.velocities_divergent];
dynamic_viscosities = zeros(numel(x_coordinates),1);
thermal_conductivities = zeros(numel(x_coordinates),1);
prandtl_numbers = zeros(numel(x_coordinates),1);
reynolds_numbers = zeros(numel(x_coordinates),1);
nusselt_numbers = zeros(numel(x_coordinates),1);
heat_transfer_coeffitients = zeros(numel(x_coordinates),1);

for i = 1:numel(x_coordinates)
    gas_results_py = combustion_module.hydro_product_properties(temperatures(i), pressures(i), composition_frac(1), composition_frac(2), composition_frac(3), composition_frac(4), composition_frac(5), composition_frac(6), composition_frac(7), composition_frac(8));
    molecular_weight_py = gas_results_py{1};
    gas_constant_py = gas_results_py{2};
    gamma_py = gas_results_py{3};
    viscosity_py = gas_results_py{4};
    thermal_conductivity_py = gas_results_py{5};
    molecular_weight_i = double(molecular_weight_py);
    gas_constant_i = double(gas_constant_py);
    gamma_i = double(gamma_py);
    dynamic_viscosities(i) = double(viscosity_py);
    thermal_conductivities(i) = double(thermal_conductivity_py);
    prandtl_numbers(i) = dynamic_viscosities(i) * (gamma_i * gas_constant_i / (gamma_i - 1)) / thermal_conductivities(i);
    reynolds_numbers(i) = densities(i) * diameters(i) * velocities(i) / dynamic_viscosities(i);
    nusselt_numbers(i) = 0.023 * (reynolds_numbers(i)^0.8) * (prandtl_numbers(i)^0.4);
    heat_transfer_coeffitients(i) = thermal_conductivities(i) * nusselt_numbers(i) / diameters(i);
end

%% Wyświetlanie wyników

disp(['Masa molowa: ', num2str(molecular_weight), ' kg/kmol']);
disp(['Stała gazowa: ', num2str(gas_constant), ' J/(kg*K)']);
disp(['Wykładnik adiabaty: ', num2str(gamma)]);
disp(['Lepkość dynamiczna: ', num2str(viscosity), ' Pa*s']);
disp(['Przewodność cieplna: ', num2str(thermal_conductivity), ' W/(m*K)']);

nexttile;
plot(x_coordinates, [Flow_parameters.mach_numbers_convergent, Flow_parameters.mach_numbers_divergent]);
xlabel("distance from entry [m]");
ylabel("mach number");
grid minor;
title("Mach Number");

nexttile;
plot(x_coordinates, pressures);
xlabel("distance from entry [m]");
ylabel("pressure [Pa]");
grid minor;
title("Pressure");

nexttile;
plot(x_coordinates, temperatures);
xlabel("distance from entry [m]");
ylabel("temperature [K]");
grid minor;
title("Temperature");

nexttile;
plot(x_coordinates, velocities);
xlabel("distance from entry [m]");
ylabel("velocity [m/s]");
grid minor;
title("Velocity");
hold off;

nexttile;
plot(x_coordinates, prandtl_numbers);
xlabel("distance from entry [m]");
ylabel("prandtl number [-]");
grid minor;
title("Prandtl number");

nexttile;
plot(x_coordinates, reynolds_numbers);
xlabel("distance from entry [m]");
ylabel("reynolds number [-]");
grid minor;
title("Reynolds number");

nexttile;
plot(x_coordinates, nusselt_numbers);
xlabel("distance from entry [m]");
ylabel("nusselt number [-]");
grid minor;
title("Nusselt number");

nexttile;
plot(x_coordinates, heat_transfer_coeffitients);
xlabel("distance from entry [m]");
ylabel("heat transfer coeffitient [W/(m^2*K)]");
grid minor;
title("Heat transfer coeffitient");

