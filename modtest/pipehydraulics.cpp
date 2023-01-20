#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <math.h>


const double GRAVITY = 9.81; // gravitational acceleration in m/s^2
const double PI = 3.14159265358979323846;
const double R = 8.3144598;
const double M_PI = 3.14159265358979323846;


// Calculates various oil pipe hydraulics using the Miller equation and MIT equation.
// The function takes the following inputs:
//   d: diameter of the pipe in inches
//   L: length of the pipe in feet
//   rho: density of the oil in lb/ft^3
//   mu: viscosity of the oil in lb*s/ft^2
//   Q: flow rate of the oil in ft^3/s
// The function returns a struct containing the following outputs:
//   Re: Reynolds number
//   f: friction factor
//   h_l: head loss in feet
//   v: velocity in ft/s
//   gpm: flow rate in gallons per minute
struct OilPipeHydraulics {
  double Re;
  double f;
  double h_l;
  double v;
  double gpm;
};

OilPipeHydraulics calculate_oil_pipe_hydraulics(double d, double L, double rho, double mu, double Q) {
  // Convert the diameter to feet
  d /= 12;

  // Calculate the Reynolds number
  double Re = rho * d * Q / mu;

  // Calculate the friction factor using the Miller equation for laminar flow (Re <= 2100)
  // and the MIT equation for turbulent flow (Re > 2100)
  double f;
  if (Re <= 2100) {
    f = 64 / Re;
  } else {
    f = 0.3164 / std::pow(std::log10(Re), 2.7);
  }

  // Calculate the head loss and velocity
  double h_l = f * L * std::pow(Q, 2) / (2 * rho * d * d * d * d);
  double v = Q / (M_PI * std::pow(d, 2) / 4);

  // Calculate the flow rate in gallons per minute
  double gpm = Q * rho / (7.48 * 8.34);

  // Return the results as a struct
  OilPipeHydraulics results;
  results.Re = Re;
  results.f = f;
  results.h_l = h_l;
  results.v = v;
  results.gpm = gpm;
  return results;
}

// Calculates various gas pipe hydraulics using flowrate models.
// The function takes the following inputs:
//   d: diameter of the pipe in inches
//   P_1: upstream pressure in psia
//   P_2: downstream pressure in psia
//   T: temperature in degrees Rankine
//   r: specific gas constant in ft*lbf/lbmol*R
//   Z: compressibility factor
//   SG: specific gravity
// The function returns a struct containing the following outputs:
//   Q: flow rate in ft^3/s
//   v: velocity in ft/s
//   h_l: head loss in ft
//   gpm: flow rate in gallons per minute
struct GasPipeHydraulics {
  double Q;
  double v;
  double h_l;
  double gpm;
};

GasPipeHydraulics calculate_gas_pipe_hydraulics(double d, double P_1, double P_2, double T, double r, double Z, double SG) {
  // Convert the diameter to feet
  d /= 12;

  // Calculate the specific weight of the gas
  double rho = P_1 * 144 / (r * T) * SG;

  // Calculate the flow rate using the Weymouth equation
  double Q_weymouth = M_PI * std::pow(d, 2) / 4 * std::sqrt(2 * r * T / (SG * Z * (P_1 + P_2)) * std::log((P_1 + P_2) / P_1));

  // Calculate the flow rate using the Modified Weymouth equation
  double Q_modified_weymouth = M_PI * std::pow(d, 2) / 4 * std::sqrt(2 * r * T / (SG * Z) * std::log((P_1 + P_2) / P_1));

  // Calculate the flow rate using the Panhandle equation
  double Q_panhandle = M_PI * std::pow(d, 2) / 4 * std::sqrt(2 * r * T / (SG * Z) * std::log(P_1 / P_2));

  // Use the flow rate model that gives the highest flow rate
  double Q = std::max({Q_weymouth, Q_modified_weymouth, Q_panhandle});

  // Calculate the velocity and head loss
  double v = Q / (M_PI * std::pow(d, 2) / 4);
  double h_l = (P_1 - P_2) / rho / 144;

  // Calculate the flow rate in gallons per minute
  double gpm = Q / 7.48 / 60;

  // Return the results as a struct
  GasPipeHydraulics results;
  results.Q = Q;
  results.v = v;
  results.h_l = h_l;
  results.gpm = gpm;
  return results;
}


// Calculates the head loss in a pipe with diameter d, length L, and roughness e.
// The flow rate is given in m^3/s, and the fluid properties (density and viscosity) are given in kg/m^3 and Pa*s, respectively.
double head_loss(double d, double L, double e, double Q, double rho, double mu) {
  double Re = 4 * Q / (PI * d * mu); // Reynolds number
  double f = 64 / Re; // Fanning friction factor
  return f * L / d * (rho * GRAVITY / 2) * std::pow(Q / (PI * d * d / 4), 2);
}


// Calculates the hydraulic pressure required to transport fluid through a pipeline.
// The function takes the following inputs:
//   d: diameter of the pipe in inches
//   L: length of the pipe in feet
//   f: Darcy friction factor
//   rho: density of the fluid in lbm/ft^3
//   g: acceleration due to gravity in ft/s^2
// The function returns the pressure in psig.
double calculate_pressure(double d, double L, double f, double rho, double g) {
  // Convert the diameter to feet
  d /=  12;

  // Calculate the pressure drop in psi
  double dP = f * L / d * rho * std::pow(2 / d, 2) / 2 / g;

  // Return the pressure in psig
  return dP + 14.7;
}


// Calculates the equivalent length and diameter ratio for a pipe fitting with loss coefficient K.
void equivalent_length_and_diameter_ratio(double K, double &Le, double &D_ratio) {
  Le = K * D_ratio;
  D_ratio = std::pow(1 / K, 0.5);
}

// Calculates the vapour pressure of a fluid with temperature T in degrees Celsius.
double vapour_pressure(double T) {
  return 6.11 * std::pow(10, 7.5 * T / (237.7 + T));
}

// Calculates the bulk modulus of a fluid with density rho and speed of sound c.
double bulk_modulus(double rho, double c) {
  return rho * c * c;
}


// Calculates the gravity effects on oil and gas pipeline hydraulics.
// The function takes the following inputs:
//   d: diameter of the pipe in inches
//   L: length of the pipe in feet
//   rho: density of the fluid in lbm/ft^3
//   SG: specific gravity of the fluid
//   g: acceleration due to gravity in ft/s^2
// The function returns the gravity head in feet.
double calculate_gravity_head(double d, double L, double rho, double SG, double g) {
  // Convert the diameter to feet
  d /= 12;

  // Calculate the gravity head
  double h_g = rho * SG * g * L / (144 * rho / SG);

  // Return the gravity head in feet
  return h_g;
}

//Calculate Pipeline Pressure Drop
double pipelinePressureDrop(double flowRate, double pipeDiameter, double length, double roughness, double density, double viscosity) {
    double reynoldsNumber = flowRate * pipeDiameter / viscosity;
    double frictionFactor = 0.0;
    if (reynoldsNumber < 2300) {
        frictionFactor = 64.0 / reynoldsNumber;
    } else {
        double d = pipeDiameter / 4.0;
        double a = log(2.51 / (reynoldsNumber * sqrt(roughness / d)));
        frictionFactor = 1.0 / pow(a, 2);
    }
    double pressureDrop = frictionFactor * length * pow(flowRate, 2) / (2 * density * pow(pipeDiameter, 5));
    return pressureDrop;
}

//Calculate Geothermal Effects
double geothermalEffect(double flowRate, double pipeDiameter, double length, double roughness, double density, double viscosity, double inletTemperature, double ambientTemperature) {
    double reynoldsNumber = flowRate * pipeDiameter / viscosity;
    double frictionFactor = 0.0;
    if (reynoldsNumber < 2300) {
        frictionFactor = 64.0 / reynoldsNumber;
    } else {
        double d = pipeDiameter / 4.0;
        double a = log(2.51 / (reynoldsNumber * sqrt(roughness / d)));
        frictionFactor = 1.0 / pow(a, 2);
    }
    double pressureDrop = frictionFactor * length * pow(flowRate, 2) / (2 * density * pow(pipeDiameter, 5));
    double heatLoss = (flowRate * (inletTemperature - ambientTemperature) * 4.186) / (density * pipeDiameter);
    return pressureDrop + heatLoss;
}

// Calculates Joule Thompsons effect (Assumes fluid to be an ideal gas)
double jouleThomsonEffect(double flowRate, double pipeDiameter, double length, double roughness, double density, double viscosity, double inletPressure, double inletTemperature, double outletPressure) {
    double reynoldsNumber = flowRate * pipeDiameter / viscosity;
    double frictionFactor = 0.0;
    if (reynoldsNumber < 2300) {
        frictionFactor = 64.0 / reynoldsNumber;
    } else {
        double d = pipeDiameter / 4.0;
        double a = log(2.51 / (reynoldsNumber * sqrt(roughness / d)));
        frictionFactor = 1.0 / pow(a, 2);
    }
    double pressureDrop = frictionFactor * length * pow(flowRate, 2) / (2 * density * pow(pipeDiameter, 5));
    double jouleThomsonCoefficient = (1 / inletTemperature) * (1 / (inletPressure - outletPressure));
    double temperatureDrop = jouleThomsonCoefficient * pressureDrop;
    return temperatureDrop;
}

//Calculating Pipe Sizing
double pipeSizing(double flowRate, double roughness, double density, double viscosity, double headLoss, double &pipeDiameter, double &length) {
    double reynoldsNumber = flowRate * pipeDiameter / viscosity;
    double frictionFactor = 0.0;
    if (reynoldsNumber < 2300) {
        frictionFactor = 64.0 / reynoldsNumber;
    } else {
        double d = pipeDiameter / 4.0;
        double a = log(2.51 / (reynoldsNumber * sqrt(roughness / d)));
        frictionFactor = 1.0 / pow(a, 2);
    }
    double velocity = (flowRate / (density * M_PI / 4 * pow(pipeDiameter, 2)));
    length = (headLoss * (2 * 9.81 * pipeDiameter)) / (frictionFactor * pow(velocity, 2));
    pipeDiameter = pow((flowRate / (velocity * M_PI / 4)), 0.5);
    return 0;
}

//Calculates Pipeline Head Curves
std::vector<double> pipelineHeadCurve(double flowRate, double pipeDiameter, double length, double roughness, double density, double viscosity, int numberOfPoints) {
    std::vector<double> headCurve;
    for (int i = 0; i <= numberOfPoints; i++) {
        double reynoldsNumber = flowRate * pipeDiameter / viscosity;
        double frictionFactor = 0.0;
        if (reynoldsNumber < 2300) {
            frictionFactor = 64.0 / reynoldsNumber;
        } else {
            double d = pipeDiameter / 4.0;
            double a = log(2.51 / (reynoldsNumber * sqrt(roughness / d)));
            frictionFactor = 1.0 / pow(a, 2);
        }
        double velocity = (flowRate / (density * M_PI / 4 * pow(pipeDiameter, 2)));
        double headLoss = (frictionFactor * pow(velocity, 2) * length) / (2 * 9.81 * pipeDiameter);
        headCurve.push_back(headLoss);
        flowRate += (flowRate / numberOfPoints);
    }
    return headCurve;
}

//Calculates for Laminar and Turbulent flow
std::string flowType(double flowRate, double pipeDiameter, double viscosity) {
    double reynoldsNumber = flowRate * pipeDiameter / viscosity;
    if (reynoldsNumber < 2300) {
        return "Laminar";
    } else {
        return "Turbulent";
    }
}

// Properties of Flowing Liquids

// Gravity and Molecular Weight
// Liquid unit measurement in kg/m^3
void fluidGravityMoleWeight(double density, double molecular_mass, double& specific_gravity, double& molecular_weight, double water_density) {
    specific_gravity = density / water_density;
    molecular_weight = molecular_mass / (1.0e-3);
}

//calculates the bubble point pressure of a liquid using the Peng-Robinson equation of state
double calcBubblePointPressure(double temperature, double acentricFactor, double criticalPressure, double criticalTemperature, double P, double T) {
    double Tr = temperature / criticalTemperature;
    double omega = acentricFactor;
    double a = (0.45724 * pow(R*criticalTemperature, 2) / criticalPressure) * (1 + (0.37464 + 1.54226 * omega - 0.26992 * pow(omega, 2)) * (1 - sqrt(Tr)));
    double b = 0.0778 * R * criticalTemperature / criticalPressure;
    double A = a * P / pow(R*T, 2);
    double B = b * P / (R*T);
    double Z = 1 + B - sqrt(pow(B, 2) - 4 * A);
    double bubblePointPressure = P * Z;
    return bubblePointPressure;
}

// Calculates Bubble Point oil formation volume factor
double calcBubblePointFVF(double pressure, double temperature, double Rs) {
    double A = 1.092 + (5.047E-5 * Rs) - (5.75E-9 * pow(Rs, 2));
    double B = 5.2E-5 + (2.4E-9 * Rs);
    double FVF = A + (B * (pressure - 14.7));
    return FVF;
}

// Calculates for Isothermal compressibility
double calcIsothermalCompressibility(double pressure, double Rs) {
    double A = 1.6E-5 + (9.2E-9 * Rs);
    double B = -1.2E-11 * Rs;
    double isothermal_compressibility = A + (B * pressure);
    return isothermal_compressibility;
}

// Calculates for Undersatured oil formation volume factor
double calcUndersaturatedFVF(double pressure, double Rs, double Rsi) {
    double A = 1 + (Rs * (1.25E-5 + (5E-9 * Rsi)));
    double B = 5.9E-5 * Rs;
    double C = -1.2E-9 * pow(Rs, 2);
    double FVF = A + (B * pressure) + (C * pow(pressure, 2));
    return FVF;
}

// Calculates for Oil density
double calcOilDensity(double api_gravity) {
    double density = 141.5 / (api_gravity + 131.5);
    return density;
}

// Calculates Dead oil viscosity
double calcDeadOilViscosity(double temperature, double pressure, double Rsi) {
    double A = pow(10, (1.7667 + 0.0601*(Rsi) - 0.0261*(pow(Rsi, 2))));
    double B = pow(10, (-11.513 - 0.744*(Rsi) + 0.139*(pow(Rsi, 2))));
    double viscosity = A * exp(B * pressure);
    return viscosity;
}

//Calculates Bubble Point Oil Viscosity
double calcBubblePointViscosity(double pressure, double temperature, double Rs) {
    double A = 1.4E-5 + (1.0E-9 * Rs);
    double B = -9.0E-13 * Rs;
    double viscosity = A + (B * pressure);
    return viscosity;
}

// Calculates for Undersaturated Oil Viscosity
double calcUndersaturatedViscosity(double pressure, double temperature, double Rs, double Rsi) {
    double A = 1 + (Rs * (1.25E-5 + (5E-9 * Rsi)));
    double B = 5.9E-5 * Rs;
    double C = -1.2E-9 * pow(Rs, 2);
    double viscosity = A + (B * pressure) + (C * pow(pressure, 2));
    return viscosity;
}

// Calculates Gas/Oil interfacial tension
double calcGasOilIFT(double oil_density, double gas_density, double oil_viscosity, double gas_viscosity, double temperature) {
    double J = (oil_density - gas_density) / (oil_density + gas_density);
    double IFT = (3.14 * oil_viscosity * gas_viscosity) / (2 * sqrt(oil_viscosity * gas_viscosity) * log(1 + (J * J)));
    return IFT;
}

// Calculates Water/oil interfacial tension
double calcWaterOilIFT(double surface_tension_oil, double surface_tension_water, double density_oil, double density_water) {
    double IFT = (surface_tension_oil - surface_tension_water) / (density_oil - density_water);
    return IFT;
}

//  calculates the size of a centrifugal pump for a pipeline system
double calcPumpSize(double flow_rate, double head_required, double specific_gravity) {
    double BEP_flow = flow_rate / pow((specific_gravity * head_required), 0.75);
    double pump_size = pow(BEP_flow, 1.0/3.0) * pow((specific_gravity * head_required), 0.25);
    // std::cout << "Pump size: " << pump_size << " cubic meters per second" << std::endl;
    return pump_size;
}

// Calculates for Heater
double calcHeaterSize(double flow_rate, double inlet_temp, double outlet_temp, double specific_heat) {
    double delta_temp = outlet_temp - inlet_temp;
    double heat_load = flow_rate * specific_heat * delta_temp;
    // std::cout << "Heater size: " << heat_load << " watts" << std::endl;
    return heat_load;
}

// Calculates for Coolers
double calcCoolerSize(double flow_rate, double inlet_temp, double outlet_temp, double specific_heat) {
    double delta_temp = inlet_temp - outlet_temp;
    double heat_removal = flow_rate * specific_heat * delta_temp;
    // std::cout << "Cooler size: " << heat_removal << " watts" << std::endl;
    return heat_removal;
}

// Calculates for heat exchangers
double calcHeatExchangerSize(double flow_rate, double inlet_temp, double outlet_temp, double specific_heat) {
    double delta_temp = abs(inlet_temp - outlet_temp);
    double heat_load = flow_rate * specific_heat * delta_temp;
    // cout << "Heat exchanger size: " << heat_load << " watts" << endl;
    return heat_load;
}

// Calculates for Mixers
double calcMixerEfficiency(double flow_rate, double mixer_size) {
    double mixer_efficiency = flow_rate / mixer_size;
    // std::cout << "Mixer size: " << mixer_size << " cubic meters per second" << std::endl;
    return mixer_efficiency;
}

// Calculates pipe couplings
double calcCouplingSize(double pipe_diameter, double pipe_wall_thickness) {
    double coupling_size = pipe_diameter + 2 * pipe_wall_thickness;
    // std::cout << "Coupling size: " << coupling_size << " meters" << std::endl;
    return coupling_size;
}

// Properties of Pipeline Fluids

//Calculation for density
double calcFluidDensity(double mass, double volume) {
    double density = mass / volume;
    return density;
}

// Calculation for viscosity
double calcFluidViscosity(double shear_stress, double shear_rate) {
    double viscosity = shear_stress / shear_rate;
    return viscosity;
}

// Calculation for fluid surface tension
double calcFluidSurfaceTension(double force, double length) {
    double surface_tension = force / length;
    return surface_tension;
}



// Calculates the kinetic energy of fluid flow in a pipe with diameter d and flow rate Q in m^3/s.
double kinetic_energy(double d, double Q) {
  return Q * Q / (2 * PI * d * d / 4);
}

// Calculates the Reynolds number for a fluid with diameter d, flow rate Q, and viscosity mu in a pipe.
double reynolds_number(double d, double Q, double mu) {
  return 4 * Q / (PI * d * mu);
}

// Calculates the friction factor for a fluid with Reynolds number Re in a pipe.
double friction_factor(double Re) {
  return 64 / Re;
}

// Calculates the temperature profile of a pipeline with length L, insulation thickness t, and outer temperature T_out.
// The inner temperature is given as a function T_in(z), where z is the distance along the pipeline.
void temperature_profile(double L, double t, double T_out, double (*T_in)(double), std::vector<double> &T) {
  T.resize(L / t + 1);
  for (int i = 0; i < T.size(); i++) {
    double z = i * t;
    T[i] = (T_out - T_in(z)) * std::exp(-z / t) + T_in(z);
  }
}

// Calculates for Pipeline Insulation
double calcInsulationThickness(double pipe_diameter, double pipe_length, double ambient_temp, double fluid_temp, double thermal_conductivity, double insulation_coefficient) {
    double Q = 2 * M_PI * pipe_length * (fluid_temp - ambient_temp) / (log(pipe_diameter / (pipe_diameter - 2 * insulation_coefficient)));
    double insulation_thickness = Q / (2 * M_PI * thermal_conductivity * (fluid_temp - ambient_temp));
    return insulation_thickness;
}

// Calculates for Fluid Compressibility
double calcFluidCompressibility(double fluid_density, double bulk_modulus) {
    double compressibility = 1 / bulk_modulus;
    return compressibility;
}

// Calculates the fluid density for a given temperature and pressure using the ideal gas law.
// The temperature is given in degrees Celsius, and the pressure is given in Pa.
double density_ideal_gas(double T, double p, double m, double R) {
return p / (R * (T + 273.15));
}



// Calculates the fluid viscosity for a given temperature using the Sutherland's law.
// The temperature is given in degrees Celsius.
double viscosity_sutherland(double T, double mu_0, double T_0, double S) {
return mu_0 * std::pow(T / T_0, 3/2) * (T_0 + S) / (T + S);
}

// Calculates the enthalpy of a fluid with temperature T in degrees Celsius and specific enthalpy h.
double enthalpy(double T, double h) {
return h + T * 4.1868;
}

// Calculates the Gibbs free energy of a fluid with temperature T in degrees Celsius, pressure p, and specific Gibbs free energy g.
double gibbs_free_energy(double T, double p, double g) {
return g + T * p / 1000;
}

// Calculates the entropy of a fluid with temperature T in degrees Celsius and specific entropy s.
double entropy(double T, double s) {
return s + T * 4.1868 / 273.15;
}

// Calculates the fugacity of a fluid with pressure p, temperature T in degrees Celsius, and specific fugacity coefficient f.
double fugacity(double p, double T, double f) {
return f * p / 1000;
}

// Calculates the partial pressure of a fluid with mole fraction x and total pressure p.
double partial_pressure(double x, double p) {
return x * p;
}

// Calculates the K-value for a fluid with temperature T in degrees Celsius and specific K-value k.
double K_value(double T, double k) {
return k * std::exp(T / 50);
}

// Calculates the Bernoulli model for a fluid with velocity v and specific Bernoulli coefficient b.
double bernoulli_model(double v, double b) {
return b * v * v / 2;
}


// Calculates the conservation of mass for a fluid with density rho, flow rate Q, and area A.
double conservation_of_mass(double rho, double Q, double A) {
return rho * Q / A;
}

// Calculates the conservation of momentum for a fluid with force F and mass flow rate rho*Q.
double conservation_of_momentum(double F, double rho, double Q) {
return F - rho * Q;
}

// Calculates Equation of state for an ideal gas
double calcIdealGasEquation(double pressure, double temperature, double molar_mass, double molar_volume) {
    double universalGasConstant = 8.314; // in J/mol·K
    double idealGasEquation = pressure * molar_volume - (universalGasConstant * temperature) / molar_mass;
    return idealGasEquation;
}


// Calculates the pressure p and temperature T in degrees Celsius for a fluid using the Peng-Robinson's equation of state.
void peng_robinson(double rho, double T, double p, double &new_rho, double &new_T, double &new_p) {
double T_c = T - 273.15;
double a = 0.45724 * std::pow(R * T_c / p, 2);
double b = 0.07780 * R * T_c / p;
double A = a * p / (std::pow(R, 2) * std::pow(T_c, 2));
double B = b * p / (R * T_c);
double alpha = (1 + std::sqrt(2)) * (1 + std::sqrt(2)) / 8;
double a_hat = alpha * A;
double b_hat = alpha * B;
double beta = b_hat / (std::pow(rho, 2) * std::pow(b, 2));
double E = -(1 - std::sqrt(2)) / 2;
double Z_c = std::pow(1 - B, 2) - (std::pow(A - 3 * std::pow(B, 2), 2) / (std::pow(A, 2)));
double Z_r = std::pow(1 - beta, 2) - (std::pow(a_hat - 3 * std::pow(b_hat, 2), 2) / (std::pow(a_hat, 2)));
double Z = Z_c + E * Z_r;
new_p = p * Z;
new_T = T * Z;
new_rho = rho * std::pow(Z, -1);
}

// Calculates the pressure p and temperature T in degrees Celsius for a fluid using the Soave-Redlich-Kwong (SRK) equation of state.
void soave_redlich_kwong(double rho, double T, double p, double &new_rho, double &new_T, double &new_p) {
double T_c = T - 273.15;
double a = 0.42747 * std::pow(R * T_c / p, 2);
double b = 0.08667 * R * T_c / p;
double A = a * p / (std::pow(R, 2) * std::pow(T_c, 2));
double B = b * p / (R * T_c);
double alpha = (1 + std::sqrt(2)) * (1 + std::sqrt(2)) / 8;
double a_hat = alpha * A;
double b_hat = alpha * B;
double beta = b_hat / (std::pow(rho, 2) * std::pow(b, 2));
double Z = 1 + (b_hat - beta - std::sqrt(std::pow(b_hat - beta, 2) + 4 * beta)) / 2;
new_p = p * Z;
new_T = T * Z;
new_rho = rho * std::pow(Z, -1);
}
