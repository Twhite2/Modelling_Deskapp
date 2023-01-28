#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <random>

double calculateNPV(double initialCost, double annualRevenue, double discountRate, double projectLife) {
    double npv = 0;
    for (int i = 0; i < projectLife; i++) {
        npv += annualRevenue / pow(1 + discountRate, i);
    }
    npv -= initialCost;
    return npv;
}

double calculatePayoutTime(double initialCost, double annualRevenue, double discountRate, double projectLife) {
    return initialCost / (annualRevenue * (1 - (1 / (pow(1 + discountRate, projectLife)))));
}



double calculateIRR(std::vector<double> cashFlows, double initialCost) {
    double irr = 0.1;
    double error = 1;
    double NPV;
    while (error > 0.0001) {
        NPV = 0;
        for (int i = 0; i < cashFlows.size(); i++) {
            NPV += cashFlows[i] / pow(1 + irr, i);
        }
        NPV -= initialCost;
        error = NPV;
        if (NPV > 0) {
            irr += 0.1 * error;
        } else {
            irr -= 0.1 * error;
        }
    }
    return irr;
}


double calculateLevelizedCost(double initialCost, double annualRevenue, double discountRate, double projectLife) {
    return initialCost + annualRevenue / (discountRate * (1 - (1 / pow(1 + discountRate, projectLife))));
}

double calculateTransportCost(double distance, double costPerMile) {
    return distance * costPerMile;
}


double monteCarloSimulation(double initialCost, double annualRevenue, double discountRate, double projectLife, int numberOfSimulations) {
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(annualRevenue, annualRevenue*0.1);
    double npv = 0;
    for (int i = 0; i < numberOfSimulations; i++) {
        double simAnnualRevenue = distribution(generator);
        npv += calculateNPV(initialCost, simAnnualRevenue, discountRate, projectLife);
    }
    return npv/numberOfSimulations;
}
