#include <iostream>
#include <random>
#include <vector>
#include <array>
#include <cmath>

class IsingModel
{
/*
    We create a class called IsingModel.
*/
private:
    int side;
    int number_sites;
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

        // Create the initial configuration of the network
        initialConfiguration();

        // Resize the neighbors vector given the number of sites
        neighbors.resize(number_sites);

        // Compute the neighbors of all sites
        computeNeighbors();

        // Precompute the probabilities
        // This depends on T, so we have to change it if we change T
        precomputeProbabilities(J, B, T);
    }

    // Create a getter for the number of sites
    int getNumberSites()
    {
        return number_sites;
    }

    // Create a function to print the network
    void printNetwork()
    {
        std::cout << "Red de Ising (" << side << "x" << side << "):\n";
        for (int i = 0; i < side; i++) {
            for (int j = 0; j < side; j++) {
                // mmapping 2D -> 1D
                int index = i * side + j;
                std::cout << network[index] << " ";
            }
            std::cout << "\n";
        }
        }

    // Create a function to set the initial configuration
    void initialConfiguration()
    {
        // For every site we choose a random number 0 or 1.
        // If it is 0 then the site takes value of -1, else of 1.
        for (int i = 0; i < number_sites; i++)
        {
            network[i] = (dist_spin(rng) == 0) ? -1 : 1;
        }
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
        // Create the array of deltaE (OK to be local)
        std::array<std::array<double, 5>, 2> deltaE;
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
        int S_index = network[idx];

        // Sum the spins of the neighbors
        int sum_neigh = 0;
        for(int n_idx : neighbors[idx]) 
        {
            sum_neigh += network[n_idx];
        }

        // Now we have to map the case with the matrix W
        // Remember that the element that we need is W[S_new == 1][sum(S_neigh)/2 + 2]
        // The condition in the row is the same as:
        int row_idx = (S_index == 1) ? 0 : 1; 
        double prob = W[row_idx][sum_neigh/2 + 2];

        if (prob >= 1.0)
        {
            network[idx] *= -1; // Accept and flip the spin
            return true;
        }
        else
        {
            // If not, we take a random number between 0 and 1.
            // If it is < prob, we flip
            if (dist_prob(rng) < prob)
            {
                network[idx] *= -1;
                return true;
            }
        }

        // If we arrive to this point, we dont flip
        return false;
    }
};

int main()
{
    // Declare the parameters of the simulation
    int side = 3;
    double J = 1.0;
    double B = 0.0;
    // Set the temperatures range
    double T_start = 4.0;
    double T_end = 1.0;
    double dT = 0.1;

    // Create the network for the fist time to reserve memory
    IsingModel model(side, J, B, T_start);
    model.printNetwork();

    // Get the number of sites
    int number_sites = model.getNumberSites();

    // We see every temperature
    for (double T = T_start; T >= T_end; T -= dT)
    {
        std::cout << "\n Simulando T = " << T << std::endl;
        // Update the temperature
        model.setTemperature(J, B, T);
        // Reset the inital configuration
        model.initialConfiguration();
        model.printNetwork();
        
        // Run metropolis
        for (int i = 0; i < number_sites * 10; i++)
        {
            model.metropolisStep();
        }
        model.printNetwork();

        // 4. Meassure
        // ...
        
    }

    return 0;
}