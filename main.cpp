#include <iostream>
#include <random>
#include <vector>
#include <array>
#include <cmath>
#include <fstream>
#include <string>

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
    IsingModel(int side, double J, double B, double T)
    {
        this->side = side;
        this->number_sites = side * side;
        // Initalize RNG
        std::random_device rd;
        rng.seed(rd());

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
    void setTemperature(double J, double B, double T) {
        precomputeProbabilities(J, B, T);
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
        double updateE = deltaE[row_idx][sum_neigh/2 + 2];

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

int main()
{
    // Declare the parameters of the simulation
    int side;
    std::cout << "Ingrese el lado: ";
    std::cin >> side;
    double J = 1.0;
    double B = 0.0;
    // Set the temperatures range
    // List with variable resolution
    std::vector<double> T_list;
    
    // High temperature
    for (double T = 5.0; T >= 2.6; T -= 0.1) T_list.push_back(T);
    
    // Critic temperature
    for (double T = 2.58; T >= 2.0 - 0.001; T -= 0.02) T_list.push_back(T);
    
    // Low temperature
    for (double T = 1.9; T >= 1.0 - 0.001; T -= 0.1) T_list.push_back(T);

    // Create the network for the fist time to reserve memory
    IsingModel model(side, J, B, T_list[0]);
    // model.printNetwork();

    // Get the number of sites
    int number_sites = model.getNumberSites();

    // Define the number of montecarlo steps
    int MCS_termalizacion = 10000;

    // Create the name of the file to store the data for thermalization
    // std::string nombre_archivo = "data_thermalization_L=" + std::to_string(side) + ".csv";
    // Open the file to save the data
    // std::ofstream archivo_term(nombre_archivo);
    // Write it
    // archivo_term << "T,MCS,Mag_per_spin,Energy_per_spin\n";

    // Create the name of the file to store the data for the steady states
    std::string nombre_res = "results_L=" + std::to_string(side) + ".csv";
    // Open the file to save the data
    std::ofstream archivo_res(nombre_res);
    // Write it
    archivo_res << "T,Mag_avg,Mag_avg2,Mag_avg4,Mag_abs,Energy_avg,Susceptibility,Specific_heat\n";
    // We see every temperature
    for (double T : T_list)
    {
        // std::cout << "\n Simulando T = " << T << std::endl;
        // Update the temperature
        model.setTemperature(J, B, T);
        // Reset the inital configuration
        model.initialConfiguration();
        // model.printNetwork();
        
        // Run montecarlo steps until thermalization
        for (int t = 0; t < MCS_termalizacion; t++)
        {
            // Run metropolis for 1 MCS
            for (int i = 0; i < number_sites; i++)
            {
                model.metropolisStep();
            }

            // // Get the data
            // double m = model.getMagnetization();
            // double E = model.getEnergy(J, B);
            // // Get the energy for spin
            // double e = E / (double)number_sites;

            // // Save it
            // archivo_term << T << "," << t << "," << m << "," << e << "\n";
        }
        // model.printNetwork();

        // Now we continue saving data after thermalization
        int number_configurations = 3000;     // How many configurations we want
        int MCS_sample = 50;      // Time between configurations

        double sum_m = 0.0, sum_m2 = 0.0, sum_m4 = 0.0, sum_m_abs = 0.0;
        double sum_e = 0.0, sum_e2 = 0.0;

        for (int m = 0; m < number_configurations; m++)
        {
            // Run the sistems for the time of decorrelation
            for (int skip = 0; skip < MCS_sample; skip++) {
                // Run 1 MCS
                for (int i = 0; i < number_sites; i++) {
                    model.metropolisStep();
                }
            }

            // Get the instantaneous m and e
            double m_inst = model.getMagnetization(); 
            double m_abs = std::abs(m_inst);
            double e_inst = model.getEnergy() / (double)number_sites;

            // Accumulate it
            sum_m  += m_inst;
            sum_m2 += m_inst * m_inst;
            sum_m4 += m_inst * m_inst * m_inst * m_inst;
            sum_m_abs += m_abs;
            sum_e  += e_inst;
            sum_e2 += e_inst * e_inst;
        }

        // Calculate the averages
        double avg_m  = sum_m / number_configurations;
        double avg_m2 = sum_m2 / number_configurations;
        double avg_m4 = sum_m4 / number_configurations;
        double avg_mabs = sum_m_abs / number_configurations;
        double avg_e  = sum_e / number_configurations;
        double avg_e2 = sum_e2 / number_configurations;

        // Calculate susceptibility and C_V
        double chi = (number_sites / T) * (avg_m2 - (avg_mabs * avg_mabs));
        double cv  = (number_sites / (T * T)) * (avg_e2 - (avg_e * avg_e));

        // Guardamos los promedios en el archivo de resultados
        archivo_res << T << "," << avg_m << "," << avg_m2 << "," << avg_m4 << "," << avg_mabs << "," << avg_e << "," << chi << "," << cv << "\n";
        
    }
    // Close the files
    // archivo_term.close();
    archivo_res.close();

    return 0;
}