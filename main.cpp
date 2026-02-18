#include <iostream>
#include <random>
#include <vector>
#include <array>

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
    std::uniform_int_distribution<int> dist_spin;
    /*
        We will create a 1D vector for the network because it is dynamic.
        To map the network with a 2x2 array we do:
        (row, column) -> index = row * L + column
    */
    std::vector<int> network;
    // Create a vector of arrays to store the neighbors of each cell  
    std::vector<std::array<int, 4>> neighbors;
public:
    // Create a constructor
    IsingModel(int side)
    {
        this->side = side;
        this->number_sites = side * side;
        // Initalize RNG
        std::random_device rd;
        rng.seed(rd());

        // Set distribution between 0 and 1
        dist_spin = std::uniform_int_distribution<int>(0, 1);

        // Resize the network using the given number of sites
        network.resize(number_sites);

        // Create the initial configuration of the network
        initialConfiguration();

        // Resize the neighbors vector given the number of sites
        neighbors.resize(number_sites);

        // COmpute the neighbors of all sites
        computeNeighbors();
    }

    // Create a getter for the network
    int getSpin(int index) const
    {
        return network[index];
    }

    // Create a function to print the network
    void printNetwork()
    {
        std::cout << "Red de Ising (" << side << "x" << side << "):\n";
        for (int i = 0; i < side; i++) {
            for (int j = 0; j < side; j++) {
                // mmapping 2D -> 1D
                int index = i * side + j;
                std::cout << getSpin(index) << " ";
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
            std::cout << i << "\n";
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

    // Create a getter that returns all the neighbors of each site
    const std::array<int, 4>& getNeighborsOf(int index) const 
    {
        return neighbors[index];
    }
    
};

int main()
{
    int side = 3;
    IsingModel model(side);

    model.printNetwork();

    return 0;
}