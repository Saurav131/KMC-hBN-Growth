#include <stdio.h>
#include <chrono>
#include <stdlib.h>
#include <math.h>
#include <time.h> 
#include <vector>
#include <algorithm>
#include "config.h"
#include "restart.h"
#define RANDOM_MAX 2147483647.0  // 2^31 - 1



// ======================================================
// Utility time functions
// ======================================================
std::chrono::high_resolution_clock::time_point get_time_now() {
    return std::chrono::high_resolution_clock::now();
}

double get_elapsed_seconds(std::chrono::high_resolution_clock::time_point start_time) {
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    return elapsed.count();
}

// ======================================================
// Event structure and global containers
// ======================================================
struct Event {
    int type;       // atom type: 1=B, 2=N
    int from_x, from_y;
    int to_x, to_y;
    int ci, cj;     // event class indices (for diagnostics)
    double rate;
};

// ======================================================
// Simulation state
// ======================================================
struct SimulationState
{
    // Counters
    int depositions = 0;
    int desorptions = 0;
    int diffusions = 0;
    int swing_events = 0;

    // Physical simulation time
    double ETS = 0.0;

    // Next deposition time
    double deposition_time = -(DEPOSITION_INTERVAL * SEED);

    // Restart information
    long long current_step = 0;
};

std::vector<Event> all_events;
double total_rate = 0.0;
long long class_count[5][5] = {0};

// ======================================================
// Helper: Add/remove event operations
// ======================================================
void add_event(int type, int fx, int fy, int tx, int ty, int ci, int cj, double rate) {
    Event e = { type, fx, fy, tx, ty, ci, cj, rate };
    all_events.push_back(e);
    total_rate += rate;
    class_count[ci][cj]++;
}

void remove_event(size_t idx) {
    if (idx >= all_events.size()) return;
    total_rate -= all_events[idx].rate;
    class_count[all_events[idx].ci][all_events[idx].cj]--;
    all_events[idx] = all_events.back();
    all_events.pop_back();
}

// ======================================================
// Localized event cleanup (removes events near a point)
// ======================================================
void remove_neighboring_events(int x, int y, int radius = 2) {
    for (size_t i = 0; i < all_events.size();) {
        int dx = abs(all_events[i].from_x - x);
        int dy = abs(all_events[i].from_y - y);
        if (dx <= radius && dy <= radius) remove_event(i);
        else i++;
    }
}

// ======================================================
// Event selection based on cumulative rate
// ======================================================
Event select_event() {
    double r2 = (double) random() / RANDOM_MAX;
    double threshold = r2 * total_rate;
    double sum = 0.0;
    for (size_t k = 0; k < all_events.size(); ++k) {
        sum += all_events[k].rate;
        if (sum > threshold) return all_events[k];
    }
    return all_events.back(); // fallback (should not happen)
}

// ======================================================
// KMC time step (Gillespie-style)
// ======================================================
double compute_kmc_timestep(double total) {
    double r = (double) random() / RANDOM_MAX;
    if (r < 1e-6) r = 1e-6; // avoid log(0)
    double timestep = -log(r) / total;
    return timestep;
}

// ======================================================
// Energy/Rate calculation helpers
// ======================================================
double compute_arrhenius_rate(double barrier, double TEMP, double freq) {
    return freq * exp(-barrier / (TEMP * BOLTZMANN_CONSTANT));
}

// ======================================================
// Local event generation logic
// ======================================================
void generate_events_near_site(int x, int y, int lattice [ROW] [COL],
    const int NN_x[6], const int NN_y[6],
    const int NNN_x[6], const int NNN_y[6],
    const double Energy_Barriers_Boron[4][4],
    const double Energy_Barriers_Nitrogen[4][4],
                           double ATTEMPT_FREQUENCY,
                           double TEMP) {
    // scan local area radius 2
    for (int dx = -2; dx <= 2; ++dx) {
        for (int dy = -2; dy <= 2; ++dy) {
            int i = (x + dx + ROW) % ROW;
            int j = (y + dy + COL) % COL;
            int stp = lattice[i][j];
            if (stp == 0) continue; // empty site
            int swings[2] = {};
            int swing[2] = {};
            int m = 0;

            // NN bonds count
            for (int k = 0; k < 3; k++) {
                int sxs = (i + ROW + NN_x[3 * (stp - 1) + k]) % ROW;
                int sys = (j + COL + NN_y[3 * (stp - 1) + k]) % COL;
                if (lattice[sxs][sys] != 0 && lattice[sxs][sys] != stp) {
                    swings[0] = sxs;
                    swings[1] = sys;
                    m++;
                }
            }

            // monomer desorption
            if (m == 0) {
                double rate = 0.0;
                double barrier = (stp == 2) ? 0.44 : 0.78;
                rate = compute_arrhenius_rate(barrier, TEMP, ATTEMPT_FREQUENCY);
                add_event(stp, i, j, i, j, 4, 3, rate);
            }

            // Diffusion / swing events
            for (int j2 = 0; j2 < 6; j2++) {
                int sx = (i + ROW + NNN_x[j2]) % ROW;
                int sy = (j + COL + NNN_y[j2]) % COL;
                if (lattice[sx][sy] != 0) continue; // target empty
                int n = 0;
                for (int k2 = 0; k2 < 3; k2++) {
                    int tsi = (sx + ROW + NN_x[3 * (stp - 1) + k2]) % ROW;
                    int tsj = (sy + COL + NN_y[3 * (stp - 1) + k2]) % COL;
                    if (lattice[tsi][tsj] != 0 && lattice[tsi][tsj] != lattice[sx][sy]) {
                        swing[0] = tsi;
                        swing[1] = tsj;
                        n++;
                    }
                }

                // swing (4,4)
                if (m == 1 && n == 1 &&
                    swings[0] == swing[0] && swings[1] == swing[1]) {
                    double barrier = (stp == 1) ? 0.32 * CORRECTION_FACTOR : 0.39 * CORRECTION_FACTOR;
                    double rate = compute_arrhenius_rate(barrier, TEMP, ATTEMPT_FREQUENCY);
                    add_event(stp, i, j, sx, sy, 4, 4, rate);
                }

                // diffusion (m,n)
                double barrier = (stp == 1) ? Energy_Barriers_Boron[m][n] : Energy_Barriers_Nitrogen[m][n];
                double rate = compute_arrhenius_rate(barrier, TEMP, ATTEMPT_FREQUENCY);
                add_event(stp, i, j, sx, sy, m, n, rate);
            }
        }
    }
}

// ======================================================
// Random deposition function 
// ======================================================
int select_deposition_site(int type, int arr[]) {
    double r = (double) random() / RANDOM_MAX;
    int site = (int) (ROW * COL * r);
    int stp = (site / 32) % 2 + 1;
    if (stp == type) {
        arr[1] = (site / 32);
        arr[0] = site % 32;
    } else {
        arr[1] = (site / 32) + 1;
        arr[0] = site % 32;
    }
    return 0;
}


int main() {
   // --------------------------------------------
   // Random SEED
   // --------------------------------------------
   srandom(12345);
   printf("Seed = %u\n", 12345);
   for (int i = 0; i < 10; i++)
       printf("%2d: %f\n", i, (double) random() / RANDOM_MAX);

   // --------------------------------------------
   // Parameters and lattice initialization
   // --------------------------------------------
   int lattice [ROW] [COL];
   for (int i = 0; i < ROW; i++)
       for (int j = 0; j < COL; j++)
           lattice[i][j] = 0;

   // --------------------------------------------
   // Simulation bookkeeping
   // --------------------------------------------
   SimulationState sim;

   long long start_iter = 0;
   int dep_saved = 0, desorp_saved = 0;
   double ETS_saved = 0.0, deptime_saved = -48.828125;
   clock_t t = clock();
   printf("Simulation Begins\n");

   auto start = get_time_now();
   bool resumed = load_restart("classhex.csv", "outhex.txt",
                               start_iter, dep_saved, desorp_saved,
                               ETS_saved, deptime_saved, lattice);

   if (resumed) {
       printf("Continuing run; coverage=%d\n", dep_saved - desorp_saved);
       sim.depositions = dep_saved;
       sim.desorptions = desorp_saved;
       sim.ETS = ETS_saved;
       sim.deposition_time = deptime_saved;
   } else {
       printf("Starting fresh\n");
   }


   all_events.clear();
   total_rate = 0.0;
   
   for(int x=0; x < ROW; x++)
   {
       for(int y=0; y < COL; y++)
       {
           if(lattice[x][y] != 0)
           {
               generate_events_near_site(
                   x,
                   y,
                   lattice,
                   NN_x,
                   NN_y,
                   NNN_x,
                   NNN_y,
                   Energy_Barriers_Boron,
                   Energy_Barriers_Nitrogen,
                   ATTEMPT_FREQUENCY,
                   TEMP
               );
           }
       }
   }

   // --------------------------------------------
   // Main KMC loop
   // --------------------------------------------
   for (long long step = start_iter + 1; step <= MAX_ITERATIONS; ++step) {

   // ---- compute deposition allowance and dep_rate for this iteration ----
   int coverage = sim.depositions - sim.desorptions;
   bool deposition_allowed = (sim.ETS > sim.deposition_time) || (coverage <= (int)SEED);

   // if there are no events and deposition is not allowed -> stop
   if (all_events.empty() && !deposition_allowed) {
       printf("No events left at step %lld\n", step);
       break;
   }

   // determine deposition type and rate (same physics as before)
   int dep_type = (sim.depositions % 2 == 0) ? 1 : 2;
   double dep_barrier = (dep_type == 1) ? 0.014 : 0.018;
   double dep_rate = 0.0;
   if (deposition_allowed) {
       dep_rate = ATTEMPT_FREQUENCY * exp(-dep_barrier / (TEMP * BOLTZMANN_CONSTANT));
   } else {
       dep_rate = 0.0;
   }

   // effective total used for timestep and selection
   double effective_total = total_rate + dep_rate;
   if (effective_total <= 0.0) {
       printf("No events (including deposition) at iteration %lld — stopping\n", step);
       break;
   }

   // draw random number and decide whether deposition wins
   double rsel = (double) random() / RANDOM_MAX;
   double Renergy = rsel * effective_total;

   if (deposition_allowed && Renergy < dep_rate) {
       // ---------- Deposition selected ----------
       int arr[2];
       select_deposition_site(dep_type, arr);            // select deposition site by your existing function
       lattice[arr[0]][arr[1]] = dep_type;       // place atom
       sim.depositions++;                                // update deposition counter
       sim.deposition_time += DEPOSITION_INTERVAL;                // same update as your original code

      // advance time using effective_total
      double timestep = compute_kmc_timestep(effective_total);
      sim.ETS += timestep;

       // Local updates around the newly deposited atom
       remove_neighboring_events(arr[0], arr[1], 2);
       generate_events_near_site(arr[0], arr[1], lattice,
                             NN_x, NN_y, NNN_x, NNN_y,
                             Energy_Barriers_Boron, Energy_Barriers_Nitrogen,
                             ATTEMPT_FREQUENCY, TEMP);

       // go to next iteration
       continue;
   }
       
       // Select & execute event
       Event ev = select_event();
       double timestep = compute_kmc_timestep(total_rate);
       sim.ETS += timestep;

       // Apply event to lattice
       lattice[ev.from_x][ev.from_y] = 0;
       lattice[ev.to_x][ev.to_y] = ev.type;

       // Update counters by class
       if (ev.ci == 4 && ev.cj == 3) { sim.desorptions++; sim.deposition_time -= DEPOSITION_INTERVAL; }
       else if (ev.ci == 4 && ev.cj == 4) { sim.swing_events++; }
       else if (ev.ci == 4 && (ev.cj == 1 || ev.cj == 2)) {
           sim.depositions++; sim.deposition_time += DEPOSITION_INTERVAL;
       } else sim.diffusions++;

       // Localized event update
       remove_neighboring_events(ev.from_x, ev.from_y, 2);
       remove_neighboring_events(ev.to_x, ev.to_y, 2);
       generate_events_near_site(ev.from_x, ev.from_y, lattice,
                             NN_x, NN_y, NNN_x, NNN_y,
                             Energy_Barriers_Boron, Energy_Barriers_Nitrogen,
                             ATTEMPT_FREQUENCY, TEMP);
       generate_events_near_site(ev.to_x, ev.to_y, lattice,
                             NN_x, NN_y, NNN_x, NNN_y,
                             Energy_Barriers_Boron, Energy_Barriers_Nitrogen,
                             ATTEMPT_FREQUENCY, TEMP);

       



        if (step % PRINT_ITERATION == 0) {

            // Print progress periodically
            double class_rate[5][5] = {0.0};
            int class_events[5][5] = {0};
            
            for (const auto &e : all_events) {
                class_rate[e.ci][e.cj] += e.rate;
                class_events[e.ci][e.cj]++;
            }


            double seconds = get_elapsed_seconds(start);
        
            FILE *fp = fopen("classhex.csv",
                             step == PRINT_ITERATION ? "w" : "a");
        
            if (fp == NULL) {
                perror("classhex.csv");
                return -1;
            }

        fseek(fp, 0, SEEK_END);

        if (ftell(fp) == 0) {
            fprintf(fp,
            "Iteration, Executed_Event_Class_i, Executed_Event_Class_j, Coverage, Total_Rate, Time_Step, Simulation_Time, Deposition_Time, Wall_Time");
        
            for (int ii = 0; ii < 5; ii++) {
                for (int jj = 0; jj < 5; jj++) {
                    if (!((ii == 0 && jj == 3) || (ii == 3 && jj == 0))) {
                        fprintf(fp,
                                ",Rate_%d_%d,Count_%d_%d",
                                ii, jj, ii, jj);
                    }
                }
            }
            fprintf(fp,"\n");
        }
        fprintf(fp,
            "%lld,%d,%d,%d,%.2e,%.2e,%.2e,%.2e,%.2e",
            step, ev.ci, ev.cj, sim.depositions-sim.desorptions, total_rate, timestep, sim.ETS, sim.deposition_time, seconds);
            
            for (int ii=0; ii<5; ii++) {
                for (int jj=0; jj<5; jj++) {
            
                    if (!((ii==0 && jj==3) || (ii==3 && jj==0))) {
            
                        fprintf(fp,
                                ",%g,%d",
                                class_rate[ii][jj],
                                class_events[ii][jj]);
            
                    }
                }
            }
            
            fprintf(fp,"\n");
            fclose(fp);}




       // Periodic structure dump
       if (step % PRINT_ITERATION == 0) {
           FILE *fp2 = fopen("outhex.txt", step == PRINT_ITERATION ? "w" : "a");
           fprintf(fp2,"ITEM: TIMESTEP\n%lld\n", step);
           fprintf(fp2,"ITEM: NUMBER OF ATOMS\n%d\n", ROW*COL);
           fprintf(fp2,"ITEM: BOX BOUNDS\n0 %d\n0 %d\n-0.5 0.5\n",6*ROW,4*COL);
           fprintf(fp2,"ITEM: ATOMS type x y z coverage\n");
           for (int b=0;b<COL;b++){
               int flip = (b % 2);
               for (int a=0;a<ROW;a++){
                   int t = lattice[a][b];
                   const char* atom = (t==1)?"B":(t==2)?"N":"C";
                   int z = (t>0)?3:0;
                   int X = (flip==1)?6*a:(6*a+3);
                   fprintf(fp2,"%s %d %d %d %d\n",atom,X,4*b,z,sim.depositions-sim.desorptions);
               }
           }
           fclose(fp2);
       }

       if (sim.depositions - sim.desorptions > 1700) break;
   }

// ======================================================
// Final Simulation Summary
// ======================================================
t = clock() - t;
double runtime = static_cast<double>(t) / CLOCKS_PER_SEC;

printf("\n==================================================\n");
printf("            KMC Simulation Complete\n");
printf("==================================================\n");

printf("Wall-clock runtime : %.2f s\n", runtime);
printf("Simulated time     : %.3e s\n", sim.ETS);

printf("\n---------------- Event Summary -------------------\n");
printf("Total active events : %zu\n", all_events.size());
printf("Total event rate    : %.3e s^-1\n", total_rate);

printf("\n---------------- Counters ------------------------\n");
printf("Depositions   : %d\n", sim.depositions);
printf("Desorptions   : %d\n", sim.desorptions);
printf("Diffusions    : %d\n", sim.diffusions);
printf("Swing events  : %d\n", sim.swing_events);
printf("Final coverage: %d\n", sim.depositions - sim.desorptions);

printf("\n----------- Active Events by Class ---------------\n");
for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
        printf("(%d,%d): %-8lld ", i, j, class_count[i][j]);
    }
    printf("\n");
}

printf("==================================================\n");

return 0;
}
