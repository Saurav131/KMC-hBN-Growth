#ifndef CONFIG_H
#define CONFIG_H

// Physical parameters (constants, not macros)
constexpr double TEMP = 1600.0;
constexpr double CORRECTION_FACTOR = 2.0;
constexpr double ATTEMPT_FREQUENCY = 1e12;  // 10^12
constexpr double BOLTZMANN_CONSTANT = 8.617e-5;  // eV/K
constexpr unsigned int SEED = 100;
constexpr int ROW = 32;
constexpr int COL = 64;
constexpr double DEPOSITION_INTERVAL = 0.48828125;

constexpr long long MAX_ITERATIONS = 100000LL;

// Output and iteration controls
constexpr long long PRINT_ITERATION = 1000; // adjust as needed

constexpr int NN_x[6] = { 0, 0, 1, 0, 0, -1 };
constexpr int NN_y[6] = { 1, -1, -1, -1, 1, 1 };
constexpr int NNN_x[6]  = { -1, 0, -1, +1, 0, +1 };
constexpr int NNN_y[6]  = { +2, +2, 0, 0, -2, -2 };

constexpr double Energy_Barriers_Boron[4][4] = {
    {0.014,0.014,0.014,0.000},
    {4.610,2.220,0.014,0.014},
    {8.170,4.530,2.450,2.140},
    {0.000,7.280,4.150,2.560}
};
constexpr double Energy_Barriers_Nitrogen[4][4] = {
    {0.018,0.018,0.018,0.000},
    {4.610,2.150,0.018,0.018},
    {5.410,4.120,2.810,1.860},
    {0.000,3.850,3.100,2.470}
};


#endif
