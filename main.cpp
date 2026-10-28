#include <iostream>
#include <random>
#include <vector>
#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdint>
#include <stdexcept>
#include <filesystem>
#include <iomanip>
#include <sstream>

class IsingModel
{
/*
    We create a class called IsingModel.
*/
private:
    int side;
    int number_sites;
    double J, B;
    // Random numbers generator
    std::mt19937 rng;
    // For initial spins
    std::uniform_int_distribution<int> dist_spin;
    // For index to visit
    std::uniform_int_distribution<int> dist_index;
    // For probability of accept the transition
    std::uniform_real_distribution<double> dist_prob;
    /*
        We will create a 1D vector for the network because it is dynamic.
        To map the network with a 2x2 array we do:
        (row, column) -> index = row * L + column
    */
    std::vector<int> network;
    // Create a vector of arrays to store the neighbors of each cell  
    std::vector<std::array<int, 4>> neighbors;
    /*
        Create an array to store the probabilities of transition
        Given the S_new and the sum(S_neigh), then the element of the
        matrix that we need is:
        W[S_new == 1][sum(S_neigh)/2 + 2]
    */
    std::array<std::array<double, 5>, 2> W;
    // And the same for the energy
    std::array<std::array<double, 5>, 2> deltaE;
    // Also we create a variable for the total energy and magnetization
    double energy;
    double magnetization;
public:
    // Create a constructor
    IsingModel(int side, double J, double B, double T, std::uint32_t seed)
    {
        if (side < 3) {
            throw std::invalid_argument("L debe ser mayor o igual a 3");
        }
        this->side = side;
        this->number_sites = side * side;
        this->J = J;
        this->B = B;
        // Reproducible random seed; pass a different seed for independent runs.
        rng.seed(seed);

        // Set distribution for spins between 0 and 1 (only int)
        dist_spin = std::uniform_int_distribution<int>(0, 1);
        // set distribution for idx between 0 and N-1
        dist_index = std::uniform_int_distribution<int>(0, number_sites - 1);
        // set distribution for accept between 0 and 1 (double)
        dist_prob = std::uniform_real_distribution<double>(0.0, 1.0);

        // Resize the network using the given number of sites
        network.resize(number_sites);

        // Resize the neighbors vector given the number of sites
        neighbors.resize(number_sites);

        // Compute the neighbors of all sites
        computeNeighbors();

        // Precompute the probabilities
        // This depends on T, so we have to change it if we change T
        precomputeProbabilities(J, B, T);

        // Create the initial configuration of the network
        initialConfiguration();
    }

    // Create a function to set the initial configuration
    void initialConfiguration()
    {
        magnetization = 0;
        // For every site we choose a random number 0 or 1.
        // If it is 0 then the site takes value of -1, else of 1.
        for (int i = 0; i < number_sites; i++)
        {
            network[i] = (dist_spin(rng) == 0) ? -1 : 1;
            // Update the magnetization
            magnetization += network[i];
        }
        // Calculate the energy
        this->energy = calculateFullEnergy(J, B);
    }

    // Ordered initial condition for thermalization tests.
    void orderedConfiguration()
    {
        for (int& s : network) s = 1;
        magnetization = number_sites;
        energy = calculateFullEnergy(J, B);
    }

    // Create a funtion to get the neighbors of each site
    // We do this once and store it in a vector
    void computeNeighbors()
    {
        // We calculate the neighbors for each site
        for (int i = 0; i < number_sites; i++)
        {
            // We take the row and column of the actual site
            int row = i / side;
            // And the column
            int col = i % side;
            // Then we take the values of the rows and columns of the 4 neighbors
            int row_down = (row + 1) % side; // and col_down = col
            int row_up = (row - 1 + side) % side; // and col_up = col
            int col_right = (col + 1) % side; // and row_right = row
            int col_left = (col - 1 + side) % side; // and row_left = row

            // Finally we use the mapping: index = row * L + column
            neighbors[i][0] = row_down * side + col;
            neighbors[i][1] = row_up * side + col;
            neighbors[i][2] = row * side + col_right;
            neighbors[i][3] = row * side + col_left;
            // for (int j = 0; j < 4; j++)
            // {
            //     std::cout << neighbors[i][j] << " ";
            // }
            // std::cout << "\n";
            
        }
    }

    // We can calculate the possible probabilities of transition
    void precomputeProbabilities(double J, double B, double T)
    {
        if (T <= 0.0) {
            throw std::invalid_argument("La temperatura debe ser positiva");
        }
        /*
            delta_E = - 2*J*S_new*sum(S_neigh) - 2*B*S_new

            We create a matrix 2x5:
            - In the first row, we put the case S_new = -1
            - In the second, S_new = 1
            - In the 5 columns, we put the different options of the first term
            which depends on the possible values of sum(S_neigh):
                - -1-1-1-1 = -4 -> deltaE[i][0]
                - -1-1-1+1 = -2 -> deltaE[i][1]
                - -1-1+1+1 = 0 -> deltaE[i][2]
                - -1+1+1+1 = 2 -> deltaE[i][3]
                - 1+1+1+1 = 4 -> deltaE[i][4]
            Note that given the S_new and the sum(S_neigh), then the element of the
            matrix that we need is:
            deltaE[S_new == 1][sum(S_neigh)/2 + 2]

            Then, we can calculate the probabilities of transition given the deltaE:
            - If deltaE[i][j] <= 0: deltaW[i][j] = 1
            - If deltaE[i][j] > 0: deltaW[i][j] = e**(-beta * deltaE[i][j])
        */
        // We use the deltaE and W global
        // Fill it in the case S_new = -1
        deltaE[0][0] = - 8.0 * J + 2.0 * B;
        deltaE[0][1] = - 4.0 * J + 2.0 * B;
        deltaE[0][2] = 2.0 * B;
        deltaE[0][3] = 4.0 * J + 2.0 * B;
        deltaE[0][4] = 8.0 * J + 2.0 * B;
        // Fill it in the case S_new = 1
        deltaE[1][0] = 8.0 * J - 2.0 * B;
        deltaE[1][1] = 4.0 * J - 2.0 * B;
        deltaE[1][2] = - 2.0 * B;
        deltaE[1][3] = - 4.0 * J - 2.0 * B;
        deltaE[1][4] = - 8.0 * J - 2.0 * B;

        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 5; j++)
            {
                if (deltaE[i][j] <= 0)
                {
                    W[i][j] = 1;
                }
                else
                {
                    W[i][j] = exp(-(1 / T) * deltaE[i][j]); // beta = 1/(kT)= 1/T if k=1
                }
            }
            
        }
        
    }

    // Method to set and update T during the simulation
    void setTemperature(double J_new, double B_new, double T) {
        // Keep the physical parameters and the energy consistent.
        this->J = J_new;
        this->B = B_new;
        energy = calculateFullEnergy(this->J, this->B);
        precomputeProbabilities(this->J, this->B, T);
    }

    // Create a function that does a metropolis step
    bool metropolisStep() 
    {
        // We begin choosing a random site
        int idx = dist_index(rng);

        // Identify the current state of the site
        int s_old = network[idx];
        int s_new = -s_old;

        // Sum the spins of the neighbors
        int sum_neigh = 0;
        for(int n_idx : neighbors[idx]) 
        {
            sum_neigh += network[n_idx];
        }

        // Now we have to map the case with the matrix W
        // Remember that the element that we need is W[S_new == 1][sum(S_neigh)/2 + 2]
        // The condition in the row is the same as:
        int row_idx = (s_new == 1) ? 1 : 0; 
        int col_idx = sum_neigh / 2 + 2;
        double prob = W[row_idx][col_idx];

        if (dist_prob(rng) < prob)
        {
            network[idx] = s_new;
            energy += deltaE[row_idx][col_idx];
            magnetization += 2 * s_new;
            return true;
        }
        return false;
    }

    // Create a function to calculate the full energy
    double calculateFullEnergy(double J, double B)
    {
        /*
            The energy is given by the following Hamiltonian:
            H = -J sum_first_neigh(S_i * S_j) - B sum(S_i)
            The first term is associated with the energy of interaction, while the second
            with the magnetic energy.
        */
        // Define the energies
        double energy_interaction = 0.0;
        double energy_magnetic = 0.0;

        // See every site
        for (int i = 0; i < number_sites; i++)
        {
            int S_i = network[i];
            
            // Magnetic term
            energy_magnetic += S_i;

            // Interaction term
            int sum_neighbors = 0;
            for (int n_idx : neighbors[i]) {
                sum_neighbors += network[n_idx];
            }
            // Here we are counting each combination twice: from i and its neighbor
            // Hence, we are going to divide by 2 (After)
            energy_interaction += S_i * sum_neighbors;
        }

        // Divide by 2 for correction
        return -0.5 * J * energy_interaction - B * energy_magnetic;
    }

    // Create a function to get the magnetization per spin
    double getMagnetization()
    {
        return magnetization / number_sites;   
    }

    // Create a function to get the energy
    double getEnergy()
    {
        return energy;
    }

    // Numerical consistency check, deliberately not called at every spin flip.
    void verifyState()
    {
        double check_m = 0.0;
        for (int s : network) {
            check_m += s;
        }
        double check_energy = calculateFullEnergy(J, B);
        if (std::abs(check_m - magnetization) > 1e-9 ||
            std::abs(check_energy - energy) > 1e-8 * number_sites) {
            throw std::runtime_error("Error: energia o magnetizacion incremental inconsistente");
        }
    }

    // Create a getter for the number of sites
    int getNumberSites()
    {
        return number_sites;
    }

    // // Create a function to print the network
    // void printNetwork()
    // {
    //     std::cout << "Red de Ising (" << side << "x" << side << "):\n";
    //     for (int i = 0; i < side; i++) {
    //         for (int j = 0; j < side; j++) {
    //             // mmapping 2D -> 1D
    //             int index = i * side + j;
    //             std::cout << network[index] << " ";
    //         }
    //         std::cout << "\n";
    //     }
    //     }
};

int main(int argc, char* argv[])
{
    // Optional seed as first command-line argument; fixed default for reproducibility.
    std::uint32_t seed = argc > 1 ? static_cast<std::uint32_t>(std::stoul(argv[1])) : 123456789u;
    std::cout << "Semilla: " << seed << "\n";
    // Quick self-test: ./ising_sim 123 --test
    if (argc > 2 && std::string(argv[2]) == "--test") {
        IsingModel check(8, 1.0, 0.0, 2.3, seed);
        check.verifyState();
        for (int t = 0; t < 100; ++t) {
            for (int i = 0; i < 8 * 8; ++i) check.metropolisStep();
            check.verifyState();
        }
        check.setTemperature(0.8, 0.2, 2.1); // Also test changes in J and B.
        check.verifyState();
        for (int t = 0; t < 100; ++t) {
            for (int i = 0; i < 8 * 8; ++i) check.metropolisStep();
            check.verifyState();
        }
        std::cout << "Consistencia de E y M verificada.\n";
        return 0;
    }
    // Preliminary time series for thermalization and autocorrelations.
    // ./ising_sim <seed> --prelim <L> <T> <random|ordered> <MCS> <save_every>
    if (argc > 2 && std::string(argv[2]) == "--prelim") {
        if (argc != 8) {
            std::cerr << "Uso: ./ising_sim <seed> --prelim <L> <T> <random|ordered> <MCS> <save_every>\n";
            return 1;
        }
        int L = std::stoi(argv[3]);
        double T = std::stod(argv[4]);
        std::string initial = argv[5];
        int total_mcs = std::stoi(argv[6]);
        int save_every = std::stoi(argv[7]);
        if (L < 3 || !std::isfinite(T) || T <= 0.0 || total_mcs <= 0 || save_every <= 0 ||
            (initial != "random" && initial != "ordered")) {
            std::cerr << "Parametros preliminares invalidos.\n";
            return 1;
        }

        IsingModel model(L, 1.0, 0.0, T, seed);
        if (initial == "ordered") model.orderedConfiguration();
        model.verifyState();
        int N = model.getNumberSites();

        std::ostringstream filename;
        filename << "prelim_L=" << L << "_T=" << std::fixed << std::setprecision(8) << T
                 << "_init=" << initial << "_seed=" << seed << ".csv";
        std::filesystem::create_directories("preliminar_data");
        auto path = std::filesystem::path("preliminar_data") / filename.str();
        if (std::filesystem::exists(path)) {
            std::cerr << "El archivo ya existe: " << path << ". Cambia la semilla o guardalo aparte.\n";
            return 1;
        }
        std::ofstream out(path);
        if (!out) throw std::runtime_error("No se pudo abrir el CSV preliminar");
        out << std::setprecision(17);
        out << "MCS,Energy_per_spin,Mag_per_spin\n";
        auto save = [&](int t) {
            out << t << "," << model.getEnergy() / N << "," << model.getMagnetization() << "\n";
        };
        save(0);  // The initial configuration, before any Monte Carlo attempt.
        for (int t = 1; t <= total_mcs; ++t) {
            for (int i = 0; i < N; ++i) model.metropolisStep();
            if (t % save_every == 0 || t == total_mcs) save(t);
        }
        model.verifyState();
        std::cout << "Datos preliminares guardados en: " << path << "\n";
        return 0;
    }

    // Production run: one (L, T, B, seed) per invocation, random initial state.
    // Usage:
    // ./ising_sim <seed> --production <L> <T> <B> <thermalization_MCS>
    //              <measurement_MCS> <sample_every_MCS> <block_MCS>
    if (argc > 2 && std::string(argv[2]) == "--production") {
        if (argc != 10) {
            std::cerr << "Uso: ./ising_sim <seed> --production <L> <T> <B> "
                         "<thermalization_MCS> <measurement_MCS> <sample_every_MCS> <block_MCS>\n";
            return 1;
        }

        int L = std::stoi(argv[3]);
        double T = std::stod(argv[4]);
        double B = std::stod(argv[5]);
        std::int64_t thermalization_mcs = std::stoll(argv[6]);
        std::int64_t measurement_mcs = std::stoll(argv[7]);
        std::int64_t sample_every = std::stoll(argv[8]);
        std::int64_t block_mcs = std::stoll(argv[9]);
        const double J = 1.0;

        if (L < 3 || !std::isfinite(T) || T <= 0.0 || !std::isfinite(B) ||
            thermalization_mcs < 0 || measurement_mcs <= 0 ||
            sample_every <= 0 || block_mcs < sample_every ||
            block_mcs % sample_every != 0 || measurement_mcs % block_mcs != 0 ||
            measurement_mcs / block_mcs < 4) {
            std::cerr << "Parametros invalidos: se requieren >=4 bloques completos; "
                         "block_MCS debe ser multiplo de sample_every_MCS y dividir measurement_MCS.\n";
            return 1;
        }

        std::ostringstream filename;
        filename << "blocks_L=" << L << "_T=" << std::fixed << std::setprecision(8) << T
                 << "_B=" << B << "_seed=" << seed << ".csv";
        std::filesystem::create_directories("data/production");
        auto path = std::filesystem::path("data/production") / filename.str();
        if (std::filesystem::exists(path)) {
            std::cerr << "El archivo ya existe: " << path << ". Usa otra semilla o renombra el archivo.\n";
            return 1;
        }

        IsingModel model(L, J, B, T, seed);
        model.verifyState();
        const int N = model.getNumberSites();
        for (std::int64_t t = 0; t < thermalization_mcs; ++t) {
            for (int i = 0; i < N; ++i) model.metropolisStep();
        }
        model.verifyState();
        std::cout << "Termalizacion completada: " << thermalization_mcs << " MCS\n";

        std::ofstream out(path);
        if (!out) throw std::runtime_error("No se pudo abrir el CSV de produccion");
        out << std::setprecision(17);
        out << "block,L,T,J,B,seed,thermalization_mcs,measurement_mcs,"
               "sample_every_mcs,block_mcs,samples,m,m_abs,m2,m4,e,e2\n";

        // Means of the sufficient moments in each time block.
        long double sum_m = 0.0L, sum_abs_m = 0.0L, sum_m2 = 0.0L, sum_m4 = 0.0L;
        long double sum_e = 0.0L, sum_e2 = 0.0L;
        std::int64_t samples = 0;
        const std::int64_t n_blocks = measurement_mcs / block_mcs;

        for (std::int64_t t = 1; t <= measurement_mcs; ++t) {
            for (int i = 0; i < N; ++i) model.metropolisStep();

            if (t % sample_every == 0) {
                double m = model.getMagnetization();
                double e = model.getEnergy() / N;
                sum_m += m;
                sum_abs_m += std::abs(m);
                sum_m2 += m * m;
                sum_m4 += (m * m) * (m * m);
                sum_e += e;
                sum_e2 += e * e;
                ++samples;
            }

            if (t % block_mcs == 0) {
                const std::int64_t block = t / block_mcs;
                if (samples != block_mcs / sample_every)
                    throw std::runtime_error("Cantidad inesperada de mediciones por bloque");

                // Each record contains means WITHIN the block, not raw sums.
                out << block << "," << L << "," << T << "," << J << "," << B << ","
                    << seed << "," << thermalization_mcs << "," << measurement_mcs << ","
                    << sample_every << "," << block_mcs << "," << samples << ","
                    << static_cast<double>(sum_m / samples) << ","
                    << static_cast<double>(sum_abs_m / samples) << ","
                    << static_cast<double>(sum_m2 / samples) << ","
                    << static_cast<double>(sum_m4 / samples) << ","
                    << static_cast<double>(sum_e / samples) << ","
                    << static_cast<double>(sum_e2 / samples) << "\n";

                sum_m = sum_abs_m = sum_m2 = sum_m4 = 0.0L;
                sum_e = sum_e2 = 0.0L;
                samples = 0;
                if (block == 1 || block % 5 == 0 || block == n_blocks)
                    std::cout << "Bloque " << block << "/" << n_blocks << " completado\n";
            }
        }
        model.verifyState();
        out.close();
        if (!out) throw std::runtime_error("Error al guardar el CSV de produccion");
        std::cout << "Datos por bloques guardados en: " << path << "\n";
        return 0;
    }

    std::cerr << "Uso:\n"
              << "  ./ising_sim <seed> --test\n"
              << "  ./ising_sim <seed> --prelim <L> <T> <random|ordered> <MCS> <save_every>\n"
              << "  ./ising_sim <seed> --production <L> <T> <B> <thermalization_MCS> "
                 "<measurement_MCS> <sample_every_MCS> <block_MCS>\n";
    return 1;
}
