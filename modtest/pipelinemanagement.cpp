#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>

// Calculates Pipeline Risk models
enum RiskLevel {LOW, MODERATE, HIGH};

RiskLevel determineRiskLevel(double pipelineAge, double pipelinePressure, double pipelineMaterial) {
    double ageFactor = 1;
    double pressureFactor = 1;
    double materialFactor = 1;
    double riskScore;
    if (pipelineAge > 20) {
        ageFactor = 1.2;
    }
    if (pipelinePressure > 100) {
        pressureFactor = 1.5;
    }
    if (pipelineMaterial == 1) {
        materialFactor = 1.1;
    }
    riskScore = ageFactor * pressureFactor * materialFactor;
    if (riskScore < 1.5) {
        return LOW;
    } else if (riskScore >= 1.5 && riskScore < 2.5) {
        return MODERATE;
    } else {
        return HIGH;
    }
}

// Calculates Pipeline Risk analysis and management

RiskLevel determineOverallRisk(std::vector<double> riskFactors) {
    double overallRisk = 1;
    for (int i = 0; i < riskFactors.size(); i++) {
        overallRisk *= riskFactors[i];
    }
    if (overallRisk < 0.2) {
        return LOW;
    } else if (overallRisk >= 0.2 && overallRisk < 0.7) {
        return MODERATE;
    } else {
        return HIGH;
    }
}

// Calculates pipeline integrity management

enum IntegrityStatus {GOOD, FAIR, POOR};

IntegrityStatus evaluateIntegrity(std::vector<double> integrityFactors) {
    double integrityScore = 0;
    for (int i = 0; i < integrityFactors.size(); i++) {
        integrityScore += integrityFactors[i];
    }
    integrityScore /= integrityFactors.size();
    if (integrityScore > 0.8) {
        return GOOD;
    } else if (integrityScore >= 0.5 && integrityScore <= 0.8) {
        return FAIR;
    } else {
        return POOR;
    }
}