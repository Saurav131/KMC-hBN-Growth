#ifndef RESTART_H
#define RESTART_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

// bring row,col from your main code
#ifndef row
#define row 32
#endif
#ifndef col
#define col 64
#endif

// small helper
static inline std::string trim(const std::string &s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

/**
 * load_restart:
 * Reads second-last entries from classhex.csv and outhex.txt
 * Restores: start_iter, dep, desorp, ETS, deptime, occ[][]
 * Returns true if successful.
 */
inline bool load_restart(const char *classfile,
                  const char *outfile,
                  long long &start_iter_out,
                  int &dep_out,
                  int &desorp_out,
                  double &ETS_out,
                  double &deptime_out,
                  int occ_out[row][col])
{
    // initialize defaults
    start_iter_out = 0;
    dep_out = 0;
    desorp_out = 0;
    ETS_out = 0.0;
    deptime_out = -48.828125;

    // --- parse classhex.csv ---
    std::ifstream ifs_csv(classfile);
    if (!ifs_csv) return false;
    std::vector<std::string> csv_lines;
    std::string line;
    while (std::getline(ifs_csv, line)) {
        line = trim(line);
        if (!line.empty()) csv_lines.push_back(line);
    }
    ifs_csv.close();
    if (csv_lines.size() < 2) return false;

    std::string target_csv = csv_lines[csv_lines.size() - 2];
    std::vector<std::string> cols;
    {
        std::stringstream ss(target_csv);
        std::string tok;
        while (std::getline(ss, tok, ',')) cols.push_back(trim(tok));
    }
    if (cols.size() < 8) return false;
    try {
        start_iter_out = std::stoll(cols[0]);
        int coverage = std::stoi(cols[3]);
        ETS_out = std::stod(cols[6]);
        deptime_out = std::stod(cols[7]);

        // Policy: dep = coverage, desorp = 0
        dep_out = coverage;
        desorp_out = 0;
    } catch (...) {
        return false;
    }

    // --- parse outhex.txt ---
    std::ifstream ifs_out(outfile);
    if (!ifs_out) return false;
    std::vector<std::string> out_lines;
    while (std::getline(ifs_out, line)) out_lines.push_back(line);
    ifs_out.close();
    if (out_lines.empty()) return false;

    std::vector<size_t> timesteps_idx;
    for (size_t i = 0; i < out_lines.size(); ++i) {
        if (trim(out_lines[i]) == "ITEM: TIMESTEP") timesteps_idx.push_back(i);
    }
    if (timesteps_idx.size() < 2) return false;

    size_t block_start = timesteps_idx[timesteps_idx.size() - 2];
    size_t atoms_line_idx = std::string::npos;
    for (size_t k = block_start; k < out_lines.size(); ++k) {
        std::string t = trim(out_lines[k]);
        if (t.rfind("ITEM: ATOMS", 0) == 0) { atoms_line_idx = k; break; }
        if (t == "ITEM: TIMESTEP" && k != block_start) break;
    }
    if (atoms_line_idx == std::string::npos) return false;

    for (int a = 0; a < row; ++a)
        for (int b = 0; b < col; ++b)
            occ_out[a][b] = 0;

    for (size_t p = atoms_line_idx + 1; p < out_lines.size(); ++p) {
        std::string s = trim(out_lines[p]);
        if (s.empty()) continue;
        if (s.rfind("ITEM:", 0) == 0) break;

        std::stringstream ss(s);
        std::string atom;
        int ax, ay, az, coverage_field;
        if (!(ss >> atom >> ax >> ay >> az >> coverage_field)) continue;

        int a_index = -1;
        if (ax % 6 == 0) a_index = ax / 6;
        else if ((ax - 3) % 6 == 0) a_index = (ax - 3) / 6;
        int b_index = -1;
        if (ay % 4 == 0) b_index = ay / 4;

        if (a_index >= 0 && a_index < row && b_index >= 0 && b_index < col) {
            int type = 0;
            if (atom == "B") type = 1;
            else if (atom == "N") type = 2;
            occ_out[a_index][b_index] = type;
        }
    }

    return true;
}

#endif
