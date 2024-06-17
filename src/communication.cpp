#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <regex>
#include <vector>
#include <array>
#include <queue>
#include <communication.hpp>

struct int3double
{
    int id;
    double position[3];
};

class Communication{
    // Subdomains update arrays to be sent to other subdomains, kept in a communication object
    // To be used in the subdomain update function
    int no_subdomains; 
    int default_size = 1; // Default size of the communication data structure
    private std::vector<int> new_particle_nos; // Number of particles communicated across the boundaries of the subdomains during the forward communication step
    // Communication steps are:
    // 1. Communicate ghost atoms; ids and positions (new included) of ghost atoms, we update all of them to avoid checking for ids


    // Data structures for forward and backward communication
    // Forward communication entails the migration of ghost atoms
    // Backward communication entails the reception of force calculation results
    // Keep an array of arrays of pointers to communication data stored in vectors of IDPOSITION structs
    private std::vector<std::array<std::vector<int3double> , 6>> forward_communication_data;       // Data structure for forward communication : no_subdomains x 26
    private std::vector<std::array<std::vector<int3double> , 6>> backward_communication_data;

public:
    Communication(int no_subdomains) : no_subdomains(no_subdomains) {
        // Initialize the communication data structures
        forward_communication_data.resize(no_subdomains);
        backward_communication_data.resize(no_subdomains);

        for (int i = 0; i < no_subdomains; i++) {
            for (int j = 0; j < 6; j++) {
                forward_communication_data[i][j] = std::vector<int3double>;
                backward_communication_data[i][j] = std::vector<int3double>;
            }
        }
        // Initialize the new particle numbers
        new_particle_nos.resize(no_subdomains, 0);
    }

    ~Communication() {}

    // Functions
    void set_forward_communication_data(int subdomain_id, int neighbouring_subdomain_id, const std::vector<int3double>& data) {
        forward_communication_data[subdomain_id][neighbouring_subdomain_id] = data;
    }

    void set_backward_communication_data(int subdomain_id, int neighbouring_subdomain_id, const std::vector<int3double>& data) {
        backward_communication_data[subdomain_id][neighbouring_subdomain_id] = data;
    }

    std::vector<int3double> get_forward_communication_data(int subdomain_id, int neighbouring_subdomain_id) const {
        return forward_communication_data[subdomain_id][neighbouring_subdomain_id];
    }

    std::vector<int3double> get_backward_communication_data(int subdomain_id, int neighbouring_subdomain_id) const {
        return backward_communication_data[subdomain_id][neighbouring_subdomain_id];
    }

    int get_forward_communication_data_size(int subdomain_id, int neighbouring_subdomain_id) const {
        return forward_communication_data[subdomain_id][neighbouring_subdomain_id].size();
    }
    int get_backward_communication_data_size(int subdomain_id, int neighbouring_subdomain_id) const {
        return backward_communication_data[subdomain_id][neighbouring_subdomain_id].size();
    }
    
    void set_new_particle_nos(int subdomain_id, int new_particle_no) {
        new_particle_nos[subdomain_id] = new_particle_no;
    }
    int get_new_particle_nos(int subdomain_id) const {
        return new_particle_nos[subdomain_id];
    }

};

