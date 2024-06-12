import cantera as ct
def calculate_combustion_products(pressure, hydrogen_fraction, oxygen_fraction):
    if hydrogen_fraction + oxygen_fraction != 1.0:
        raise ValueError("Frakcje wodoru i tlenu muszą sumować się do 1.")
    gas = ct.Solution('gri30.yaml')
    composition = {'H2': hydrogen_fraction, 'O2': oxygen_fraction}
    gas.TPX = 600.0, pressure, composition  # Startowa temperatura 300 K, zadane ciśnienie
    reactor = ct.IdealGasConstPressureReactor(gas)
    sim = ct.ReactorNet([reactor])
    sim.advance_to_steady_state()
    gas = reactor.thermo
    temperature = gas.T
    pressure = gas.P
    composition = gas.mole_fraction_dict()
    return temperature, pressure, composition

def hydro_product_properties(temperature, pressure, c_H, c_H2, c_H2O, c_H2O2, c_HO2, c_O, c_O2, c_OH):
    composition = {'H': c_H, 'H2': c_H2, 'H2O': c_H2O, 'H2O2': c_H2O2, 'HO2': c_HO2, 'O': c_O, 'O2': c_O2, 'OH': c_OH}
    gas = ct.Solution('gri30.yaml')
    gas.TPX = temperature, pressure, composition
    molecular_weight = gas.mean_molecular_weight
    gas_constant = ct.gas_constant / molecular_weight
    gamma = gas.cp / gas.cv
    viscosity = gas.viscosity
    thermal_conductivity = gas.thermal_conductivity
    return molecular_weight, gas_constant, gamma, viscosity, thermal_conductivity