#include <cmath>
#include <math.h>
#include <iostream>

// This function calculates the liquid holdup in a pipeline using the formula:
// holdup = flow rate / (area x density x superficial velocity) where area is calculated using pipe diameter,
// density is the density of the liquid flowing through the pipe and superficial velocity is the velocity of the fluid.
double liquid_holdup(double flow_rate, double pipe_diameter, double liquid_density,double superficial_velocity) {
    double area = 3.14159 * (pipe_diameter / 2) * (pipe_diameter / 2);
    return (flow_rate / (area * liquid_density * superficial_velocity));
}

// Calculation for Phase Fraction
//total_flow_rate: gas + liquid + solid and so on
double phase_fraction(double flow_rate_phase, double total_flow_rate) {
    return (flow_rate_phase / total_flow_rate);
}

// Calculates for Slip Velocity
double slipVelocity(double density_1, double density_2, double velocity_1, double velocity_2) {
    //Calculate the average density of the two phases
    double avgDensity = (density_1 + density_2) / 2;
    //Calculate the slip velocity
    double slipVel = std::sqrt(avgDensity * (velocity_1 - velocity_2) * (velocity_1 - velocity_2) / (density_1 * density_2));
    return slipVel;
}

// Calculates for flow patterns using parameters

enum FlowPattern {SLUG, BUBBLE, STRATIFIED, CHURN, ANNULAR, DISPERSED};

FlowPattern determineFlowPattern(double liquidVelocity, double gasVelocity, double liquidFraction) {
    if (liquidVelocity > gasVelocity && liquidFraction > 0.9) {
        return SLUG;
    } else if (liquidVelocity < gasVelocity && liquidFraction < 0.1) {
        return BUBBLE;
    } else if (liquidVelocity > gasVelocity && liquidFraction > 0.1 && liquidFraction < 0.9) {
        return STRATIFIED;
    } else if (liquidVelocity > gasVelocity && liquidFraction > 0.1 && liquidFraction < 0.9) {
        return CHURN;
    } else if (liquidVelocity < gasVelocity && liquidFraction > 0.9) {
        return ANNULAR;
    } else {
        return DISPERSED;
    }
}

// Calculates for Flow Regime

enum FlowRegime {LAMINAR, TRANSITIONAL, TURBULENT};

FlowRegime determineFlowRegime(double ReynoldsNumber, double voidFraction) {
    if (ReynoldsNumber < 2300 && voidFraction < 0.3) {
        return LAMINAR;
    } else if (ReynoldsNumber >= 2300 && ReynoldsNumber < 4000 && voidFraction > 0.3 && voidFraction < 0.7) {
        return TRANSITIONAL;
    } else {
        return TURBULENT;
    }
}

// Calculates for State Flows
enum FlowState {STEADY, TRANSIENT, PSEUDO_STEADY};

FlowState determineFlowState(double timeConstant, double flowRate) {
    if (timeConstant < 1) {
        return STEADY;
    } else if (timeConstant >= 1 && timeConstant <= 10) {
        return PSEUDO_STEADY;
    } else {
        return TRANSIENT;
    }
}

