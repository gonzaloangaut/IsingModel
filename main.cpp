#include <iostream>
#include <random>
#include <vector>

class IsingModel
{
/*
    We create a class called IsingModel.
*/
private:
    int side;
    int number_sites;
    /*
        We will create a 1D vector for the network because it is dynamic.
        To map the network with a 2x2 array we do:
        (row, column) -> index = row * L + column
    */
    // std::vector<int> network;
    // Random numbers generator
    std::mt19937 rng;
    std::uniform_int_distribution<int> dist_spin;
    std::vector<int> network;  
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

        // Resize the network
        network.resize(number_sites);    
    }

    // Create a getter for the network
    int getSpin(int index) const
    {
        return network[index];
    }

    void initialConfiguration()
    {
        // For every site we choose a random number 0 or 1.
        // If it is 0 then the site takes value of -1, else of 1.
        for (int i = 0; i < number_sites; i++)
        {
            network[i] = (dist_spin(rng) == 0) ? -1 : 1;
        }
    }
};

int main()
{
    int side = 2;
    IsingModel model(side);

    model.initialConfiguration();
    for(int i = 0; i<4; i++){
        std::cout << model.getSpin(i) << "\n";
    }

    return 0;
}