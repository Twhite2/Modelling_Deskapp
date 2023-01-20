#include <iostream>
#include <cmath>
#include <math.h>


const double M_PI = 3.14159265358979323846;

// Calculation for CO2 corrosion rate

double calcCO2Corrosion(double temperature, double temperature_coefficient, double pressure, double pressure_coefficient, double CO2_mole_composition, double CO2_mole_composition_coefficient, double fugacity, double fugacity_coefficient, double pH, double pH_coefficient, double glycol_addition, double glycol_addition_coefficient) {
    // Corrosion rate calculation based on the parameters
    double corrosion_rate = temperature_coefficient * temperature + pressure_coefficient * pressure + CO2_mole_composition_coefficient * CO2_mole_composition + fugacity_coefficient * fugacity + pH_coefficient * pH + glycol_addition_coefficient * glycol_addition;
    return corrosion_rate;
}

// Calculate for H2S corrosion rate
double calcH2SCorrosion(double temperature, double temperature_coefficient, double pressure, double pressure_coefficient, double H2S_mole_composition, double H2S_mole_composition_coefficient, double fugacity, double fugacity_coefficient, double pH, double pH_coefficient, double glycol_addition, double glycol_addition_coefficient) {
    // Corrosion rate calculation based on the parameters
    double corrosion_rate = temperature_coefficient * temperature + pressure_coefficient * pressure + H2S_mole_composition_coefficient * H2S_mole_composition + fugacity_coefficient * fugacity + pH_coefficient * pH + glycol_addition_coefficient * glycol_addition;
    return corrosion_rate;
}

// Calculates for pipeine erosion (kg/m^2s)
double calcPipelineErosion(double fluid_velocity, double fluid_density) {
    // Erosion rate calculation based on the fluid velocity and density
    double erosion_rate = fluid_velocity * fluid_density;
    return erosion_rate;
}

// Get temperature critical
double getTemperatureCritical(double gas_composition) {
    // function that returns the critical temperature of the gas based on its composition
    double Tc;
    // code for calculating Tc based on the gas composition
    return Tc;
}

// Get Pressure critical
double getPressureCritical(double gas_composition) {
    // function that returns the critical pressure of the gas based on its composition
    double Pc;
    // code for calculating Pc based on the gas composition
    return Pc;
}


// Calculates for pipeline hydrates
bool calcHydrateFormation(double temperature, double pressure, double gas_composition) {
    // Hydrate formation conditions calculation based on the temperature, pressure and gas composition
    double temperature_crit = getTemperatureCritical(gas_composition); // function that returns the critical temperature of the gas
    double pressure_crit = getPressureCritical(gas_composition); // function that returns the critical pressure of the gas
    double hydrate_formation = false;
    if (temperature <= temperature_crit && pressure <= pressure_crit) {
        hydrate_formation = true;
    }
    return hydrate_formation;
}

// Calculates for Pipeline Slug Analysis
void calcPipelineSlugAnalysis(double flow_rate, double pressure, double pipeline_diameter, double liquid_density, double liquid_velocity_critical, double pressure_drop) {
    // Slug flow frequency and magnitude calculation based on the flow rate, pressure, pipeline diameter and liquid density
    double slug_frequency = 0;
    double slug_magnitude = 0;
    double flow_velocity = flow_rate / (pipeline_diameter * pipeline_diameter * M_PI / 4);
    double liquid_velocity = liquid_density * flow_velocity;
    if (liquid_velocity > liquid_velocity_critical) {
        slug_frequency = 1;
        slug_magnitude = pressure_drop / liquid_density;
    }
    std::cout << "Slug flow frequency: " << slug_frequency << std::endl;
    std::cout << "Slug flow magnitude: " << slug_magnitude << std::endl;
}

// Get Wax Appearance
double getWaxAppearanceTemperature(double oil_composition) {
    // function that returns the wax appearance temperature of the oil based on its composition
    double wax_appearance_temp;
    // code for calculating the wax appearance temperature based on the oil composition
    return wax_appearance_temp;
}

// Get wax deposite
double getWaxDepositionRate(double flow_rate, double temperature) {
    // function that returns the wax deposition rate based on the flow rate and temperature of the fluid
    double wax_deposition_rate;
    // code for calculating the wax deposition rate based on the flow rate and temperature
    return wax_deposition_rate;
}


// Calculates pipeline wax formation and deposit
bool calcWaxFormation(double temperature, double flow_rate, double oil_composition) {
    // Wax formation and deposition calculation based on the temperature, flow rate, and oil composition
    bool wax_formation = false;
    double wax_appearance_temp = getWaxAppearanceTemperature(oil_composition); // function that returns the wax appearance temperature of the oil
    double wax_deposition_rate = getWaxDepositionRate(flow_rate, temperature); // function that returns the wax deposition rate based on flow rate and temperature
    if (temperature <= wax_appearance_temp && wax_deposition_rate > 0) {
        wax_formation = true;
    }
    return wax_formation;
}

// Pipe wall shear stress
double calcPipeWallShearStress(double fluid_velocity, double fluid_viscosity, double pipe_diameter) {
    // Pipe wall shear stress calculation based on fluid velocity, fluid viscosity, and pipe diameter
    double shear_stress = (fluid_velocity * fluid_viscosity) / (2 * pipe_diameter);
    return shear_stress;
}
 