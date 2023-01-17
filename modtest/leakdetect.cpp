#include <cmath>
#include <algorithm>
const double M_PI = 3.14159265358979323846;

//Leak Detection

// Detects a leak in an oil or gas pipeline with length L and wall thickness t.
// The pipeline is assumed to be at atmospheric pressure and the leak is assumed to be an orifice flow.
// The fluid properties (density and viscosity) are given in kg/m^3 and Pa*s, respectively.
// The pipeline is monitored at regular intervals of time t_interval.
// The function returns true if a leak is detected, and false otherwise.
bool leak_detection(double L, double t, double rho, double mu, double t_interval) {
  // Set the threshold for the leak size that would trigger an alarm
  double s_threshold = 0.001; // 1 mm

  // Set the initial pressure and flow rate to be equal to the atmospheric pressure and zero, respectively
  double p_i = 101325; // atmospheric pressure in Pa
  double Q = 0;

  // Monitor the pipeline at regular intervals of time
  for (double t_elapsed = 0; t_elapsed < t_interval; t_elapsed += 0.1) {
    // Calculate the pressure and flow rate at the current time
    double p = p_i - rho * std::pow(Q / (M_PI * std::pow(t, 2) / 4), 2) / (2 * t);
    double Q_new = std::sqrt((p_i - p) * M_PI * std::pow(t, 2) * 2 / rho);

    // Check if the flow rate has increased significantly, indicating a leak
    if (std::abs(Q_new - Q) > s_threshold) {
      return true; // leak detected
    }

    // Update the pressure and flow rate for the next iteration
    p_i = p;
    Q = Q_new;
  }

  return false; // no leak detected
}


// Calculates the distance from the start of a pipeline to the location of a leak with size s and flow rate Q.
// The pipeline is assumed to be at atmospheric pressure and the leak is assumed to be an orifice flow.
// The fluid properties (density and viscosity) are given in kg/m^3 and Pa*s, respectively.
double leak_localization(double s, double Q, double rho, double mu) {
  double d = 2 * s; // diameter of the leak
  double v = Q / (M_PI * std::pow(d, 2) / 4); // velocity of the fluid at the leak
  double x = mu * s / (rho * v * d); // distance from the start of the pipeline to the location of the leak
  return x;
}

// Calculates the time it takes for a leak to be detected in a pipeline with length L, wall thickness t, and leak size s.
// The pipeline is assumed to be at atmospheric pressure and the leak is assumed to be an orifice flow.
// The fluid properties (density and viscosity) are given in kg/m^3 and Pa*s, respectively.
double time_of_leak_detection(double L, double t, double s, double rho, double mu) {
  double d = 2 * t + s; // diameter of the leak
  double Q = mu * s / (rho * d * t); // leak flow rate
  double t_detection = L / Q; // time it takes for the leak to be detected
  return t_detection;
}

// Calculates the pressure at the point of leak in a pipeline with internal pressure p_i and diameter d.
// The pipeline is assumed to be at atmospheric pressure and the leak is assumed to be an orifice flow.
// The fluid properties (density and viscosity) are given in kg/m^3 and Pa*s, respectively.
double pressure_at_leak_point(double p_i, double d, double rho, double mu) {
double A = M_PI * std::pow(d, 2) / 4; // cross-sectional area of the pipeline
double v = std::sqrt(2 * p_i / rho); // velocity of the fluid in the pipeline
double p_l = p_i + rho * v * v / 2 - rho * std::pow(4 * mu * v / (rho * A * d), 2) / (2 * p_i); // pressure at the point of leak
return p_l;
}

// Calculates the leak volume in a pipeline with length L, wall thickness t, and leak size s.
// The pipeline is assumed to be at atmospheric pressure and the leak is assumed to be an orifice flow.
// The fluid properties (density and viscosity) are given in kg/m^3 and Pa*s, respectively.
double leak_volume(double L, double t, double s, double rho, double mu) {
double d = 2 * t + s; // diameter of the leak
double Q = mu * s / (rho * d * t); // leak flow rate
double V = Q * L / t; // leak volume
return V;
}

//Calculates Gas Leak Detection

// Detects a gas leak in a pipeline with length L and wall thickness t.
// The pipeline is assumed to be at atmospheric pressure and the leak is assumed to be an orifice flow.
// The fluid properties (density and viscosity) are given in kg/m^3 and Pa*s, respectively.
// The pipeline is monitored at regular intervals of time t_interval.
// The function returns true if a leak is detected, and false otherwise.
bool gas_leak_detection(double L, double t, double rho, double mu, double t_interval) {
  // Set the threshold for the leak size that would trigger an alarm
  double s_threshold = 0.001; // 1 mm

  // Set the initial pressure and flow rate to be equal to the atmospheric pressure and zero, respectively
  double p_i = 101325; // atmospheric pressure in Pa
  double Q = 0;

  // Monitor the pipeline at regular intervals of time
  for (double t_elapsed = 0; t_elapsed < t_interval; t_elapsed += 0.1) {
    // Calculate the pressure and flow rate at the current time
    double p = p_i - rho * std::pow(Q / (M_PI * std::pow(t, 2) / 4), 2) / (2 * t);
    double Q_new = std::sqrt((p_i - p) * M_PI * std::pow(t, 2) * 2 / rho);

    // Check if the flow rate has increased significantly, indicating a leak
    if (std::abs(Q_new - Q) > s_threshold) {
      return true; // leak detected
    }

    // Update the pressure and flow rate for the next iteration
    p_i = p;
    Q = Q_new;
  }

  return false; // no leak detected
}


//Calculates Liquid Leak Detection

// Detects a liquid leak in an oil pipeline with length L and wall thickness t.
// The pipeline is assumed to be at atmospheric pressure and the leak is assumed to be an orifice flow.
// The fluid properties (density and viscosity) are given in kg/m^3 and Pa*s, respectively.
// The pipeline is monitored at regular intervals of time t_interval.
// The function returns true if a leak is detected, and false otherwise.
bool liquid_leak_detection(double L, double t, double rho, double mu, double t_interval) {
  // Set the threshold for the leak size that would trigger an alarm
  double s_threshold = 0.001; // 1 mm

  // Set the initial pressure and flow rate to be equal to the atmospheric pressure and zero, respectively
  double p_i = 101325; // atmospheric pressure in Pa
  double Q = 0;

  // Monitor the pipeline at regular intervals of time
  for (double t_elapsed = 0; t_elapsed < t_interval; t_elapsed += 0.1) {
    // Calculate the pressure and flow rate at the current time
    double p = p_i - rho * std::pow(Q / (M_PI * std::pow(t, 2) / 4), 2) / (2 * t);
    double Q_new = std::sqrt((p_i - p) * M_PI * std::pow(t, 2) * 2 / rho);

    // Check if the flow rate has increased significantly, indicating a leak
    if (std::abs(Q_new - Q) > s_threshold) {
      return true; // leak detected
    }

    // Update the pressure and flow rate for the next iteration
    p_i = p;
    Q = Q_new;
  }

  return false; // no leak detected
}


// Detect Full and Partial leak opening

// Detects whether a leak in a pipeline with size s and flow rate Q is fully or partially open.
// The pipeline is assumed to be at atmospheric pressure and the leak is assumed to be an orifice flow.
// The fluid properties (density and viscosity) are given in kg/m^3 and Pa*s, respectively.
// The function returns true if the leak is fully open, and false if the leak is partially open.
bool leak_opening_detection(double s, double Q, double rho, double mu) {
  // Calculate the flow rate through the leak at full opening
  double Q_full = rho * std::pow(s, 2) / (2 * mu);

  // Check if the flow rate is equal to the flow rate at full opening, indicating a fully open leak
  if (std::abs(Q - Q_full) < 1e-6) {
    return true; // fully open leak
  } else {
    return false; // partially open leak
  }
}


// Calculates the orifice flow coefficient for a leak with size s and diameter d.
double orifice_flow_coefficient(double s, double d) {
double C_d = s / d; // orifice flow coefficient
return C_d;
}

// Calculates the leak containment in a pipeline with length L, wall thickness t, and leak size s.
// The pipeline is assumed to be at atmospheric pressure and the leak is assumed to be an orifice flow.
// The fluid properties (density and viscosity) are given in kg/m^3 and Pa*s, respectively.
// The containment efficiency e is given as a fraction (e.g. 0.9 for 90% efficiency).
double leak_containment(double L, double t, double s, double rho, double mu, double e) {
double V = leak_volume(L, t, s, rho, mu); // leak volume
double V_containment = e * V; // leak containment volume
return V_containment;
}