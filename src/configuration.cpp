#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP   

#include "configuration.hpp"
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

// The configuration object keeping track of the current configuration
// of the program.

struct IDPOSITION
{
    int id;
    double position[3];
};


/***** Private *****/


/***** Public *****/
namespace {
    // Local functions


}



class Configuration{
    public: 
    int no_cores;
    int no_subdomains; // number of subdomains - these are high level regions, that are assigned a processor owning it, procesors then iterate over owned regions
    int no_particles; 
    double box_size;
    int no_iterations;
    double subdomain_size;
    double dt;
    int no_subdomain_particles;

    // Subdomain boundaries and neighbours
    std::array<std::array<int, this->no_subdomains>, 6> subdomain_boundaries;       
    std::array<std::array<int, 26>, this->no_subdomains> subdomain_neighbours;
    std::array<ProcessorSubdomain, this->no_subdomains^3> processor_subdomains;
    Communication communication;


       Configuration(int no_cores, 
                int no_subdomains,  
                double dt, 
                int no_iterations, 
                double box_size, 
                int no_particles,
                int no_subdivisions){
        /**
         * @param no_cores The number of cores to be used.
         * @param no_subdivisiopn The number of subdivisions.
         * @param dt The time step.
         * @param no_iterations The number of iterations.
         * @param box_size The size of the box.
         * @param no_particles The number of particles.
         */

        this->no_cores = no_cores;
        this->no_subdivisions = no_subdivisions;
        this-> dt = dt;
        this-> no_iterations = no_iterations;
        this-> box_size = box_size;
        this-> no_particles = no_particles;
        this->no_particles_subdomain = no_particles/no_subdomains;

    }




    int init(){
        // Initialize the subdomain boundaries
        init_subdomain_boundaries(); 
        // Initialize the subdomain neighbours
        init_subdomain_neighbours();
        // Parallell initialize the subdomains that are assigned to processors
        // #pragma
        for (int i = 0; i < this->no_subdomains^3; i++){
            // TODO : arguments, subdomain size
            this->processor_subdomains[i] = ProcessorSubdomain();
        }
    }
    void init_subdomain_boundaries(){
        /**
         * Initializes the subdomains that are assigned to processors.
         * The cuboid simulation region is divided into no_subdomains^3 cubes.
         * Each row of the array represents a subdomain, with bounds given by (xl, xu, yl, yu, zl, zu).
         *
         * @throws None
         */

        double xl; double yl; double zl;
        double xu; double yu; double zu;
        double sz = this->subdomain_size;

        for (int i = 0; i < this->no_subdomains; i++){
            for (int j = 0; j < this->no_subdomains; j++){
                for (int k = 0; k < this->no_subdomains; k++){
                    // Create a subdomain - the first is given by (0, A, 0, A, 0, A)
                    xl = sz*i; 
                    yl = sz*j;
                    zl = sz*k;

                    xu = xl + sz;
                    yu = yl + sz;
                    zu = zl + sz;

                    this-> subdomain_boundaries[i][0] = xl;
                    this-> subdomain_boundaries[i][1] = xu;
                    this-> subdomain_boundaries[i][2] = yl;
                    this-> subdomain_boundaries[i][3] = yu;
                    this-> subdomain_boundaries[i][4] = zl;
                    this-> subdomain_boundaries[i][5] = zu;

                }
            }
        }
    }


    void init_neighbouring_subdomains(){
        // Initialize the neighbouring subdomain list of each subdomain
        for (int i = 0; i < this->no_subdivisions; i++){
            for (int j = 0; j < this->no_subdivisions; j++){
                for (int k = 0; k < this->no_subdivisions; k++){
                    int cell_id = this->subdomain_id_to_linear_id(i, j, k);
                    this->neighbouring_subdomains[cell_id] = this->get_neighbouring_subdomains(i, j, k);
                }
            }
        }
    }

    std::array<int, 26> get_neighbouring_subdomains(int i, int j, int k){
        // Initialize the neighbour subdomain list of each subdomain
        // Take modulo of iu, and if id < 0 wrap back to no_subdivisions-1
        int iu = i+1 == this->no_subdomains ? 0 : i+1;
        int id = i-1 == -1 ? this->no_subdomains-1 : i-1;
        int ju = j+1 == this->no_subdomains ? 0 : j+1;
        int jd = j-1 == -1 ? this->no_subdomains-1 : j-1;
        int ku = k+1 == this->no_subdomains ? 0 : k+1;
        int kd = k-1 == -1 ? this->no_subdomains-1 : k-1;

        std::array<int, 26> neighbouring_subdomains;
        // Initialize all 26 combinations
        int idx = 0; 
        for (int di = -1; di <= 1; ++di) {
            for (int dj = -1; dj <= 1; ++dj) {
                for (int dk = -1; dk <= 1; ++dk) {
                    if (di == 0 && dj == 0 && dk == 0) continue; // Skip the current cell

                    // If an index is -1, set to no_subdivisions-1, if it is no_subdivisions, set to 0

                    int ni = (i + di + this->no_subdomains) % this->no_subdomains;
                    int nj = (j + dj + this->no_subdomains) % this->no_subdomains;
                    int nk = (k + dk + this->no_subdomains) % this->no_subdomains;

                    neighbouring_subdomains[idx++] = cell_id_to_linear_id(ni, nj, nk);
                }
            }
        }

        return neighbouring_cells;

    }

    void init_processor_subdomain_ownership(){
        // Initialize processor ownership subdomain ownership
        int j = 0;
        for (int i = 0; i < this->no_cores; i++){
            int k = 0;
            for (j < this->(i+1)*this->subdomains_per_core; j++, k++){
                this->processor_subdomain_ownership[i][k] = i;
            }
        }
    }

}



class Configuration{
    public: 

    // Variables
    int no_cores;
    int no_subdivisions;
    int no_subdomains; // number of subdomains - these are high level regions, that are assigned a processor owning it, procesors then iterate over owned regions
    double dt; // time step
    int no_iterations;
    double box_size;
    int dump_step;
    int no_particles;
    int iteration = 0;



    // Misc.
    double subdomain_size;
    int subdomains_per_core; // number of subdomains per core

    // Data structures required:
    // Cell lists
    // Neighbour lists



    // Subdomain boundaries and neighbours
    std::array<std::array<int, this->no_subdomains>, 6> subdomain_boundaries;
    std::array<std::array<int, 26>, this->no_subdomains> subdomain_neighbours;

    std::vector<std::pair<int *, size_t> > neighbor_lists;   // We keep track of the neighbor lists through a vector of pointers to the neighbor lists
    std::vector<std::vector<int> > processor_ids;            // Keep track of the ids of cells that are owned by each processor

    // Position, velocity and acceleration arrays
    std::array<double, 3> positions;
    std::array<double, 3> velocities;
    std::array<double, 3> accelerations;

    // Communication structures; to be implemented
    // Each processor keeps a queue of particles to integrate on its next step

    std::vector<std::queue<int> > processor_communication_queue;
    std::vector<std::queue<int> > processor_ownership_queue;
    std:array<std::array<int, this->no_cores>, this->subdomains_per_core> processor_subdomains;

    std::queue<int> processing_queue; // Queue of subdomains to process, processor synchronously read one subdomain at a time and 
                                        // the processed subdomains are removed from the queue

    // Ghost atoms, each processor keeps track of its atoms that are ghost atoms to other processors
    // Orientations are indicated by three pairs of binary numbers: 10, 10, 10, indicating translation to sides, depthwise, and heightwise


    // When updating the cell grid we can check if the particles move sufficiently far such that they can cross a cell boundary
    std::vector<bool> update_cell_grid;
    std::vector<double> cell_grid_min_distances;
    std::vector<std::pair<int *, size_t> > cell_lists;       // We keep track of the cell lists through a vector of pointers to the cell lists
    std::vector<std::array<int, 26>> neighbouring_cells;       // Each cell keeps the ids of its neighbours

    // IDEAS: Processors calc. in a checkerbord, force calculations upon ghost atoms are computed and communicated to other processors




    int run(){
        // Run the MD simulation for the specified number of iterations
        // and return the final configuration.
        for (int i = 0; i < this->no_iterations; i++){
            // Run the MD simulation for the specified number of iterations
            // and return the final configuration.

            // Update the neighbors, through first the cell list and then the neighbor list
            // Calculate the forces
            // Update the positions etc, 
            // Communicate particles
            //update_ownership(); 

            //update_cell_lists();

            //update_neighbor_lists();

            //update_forces();

            //update_positions();

            //communicate_particles();

            // Dump the configuration at the specified dump step
            //dump();
            

            this->iteration = i;

        }

    }

    int init(){

        // Particles
        init_particle_positions();
        init_particle_velocities();

        // Subdomains
        init_subdomains();
        init_neighbouring_subdomains();
        init_procesor_subdomain_ownership();


        // Cell grid
        init_cell_lists();
        init_cell_grid();

        // Neighbor lists
        init_neighbor_lists();

        // Communication structures

        // Ghost atoms

        //



    }

    std::tuple<int, int, int> linear_id_to_subdomain_id(int lid){
        /**
         * Converts a linear index to a tuple of subdomain indices.
         *
         * @param lid the linear index to convert
         *
         * @return a tuple of three integers representing the subdomain indices (i, j, k)
         *
         * @throws None
         */
        int n = this->no_subdomains;
        int i = lid / (n * n);
        int j = (lid % (n * n)) / n;
        int k = lid % n;
        return std::make_tuple(i, j, k);
    }

    std::tuple<int, int, int> linear_id_to_cell_id(int lid){
        /**
         * Converts a linear index to a tuple of cell indices.
         *
         * @param lid the linear index to convert
         *
         * @return a tuple of three integers representing the cell indices (i, j, k)
         *
         * @throws None
         */
        int n = this->no_subdivisions;
        int i = lid / (n * n);
        int j = (lid % (n * n)) / n;
        int k = lid % n;
        return std::make_tuple(i, j, k);
    }

    int cell_id_to_linear_id(int i, int j, int k){
        /**
         * Converts a tuple of cell indices to a linear index.
         *
         * @param i the first cell index
         * @param j the second cell index
         * @param k the third cell index
         *
         * @return the linear index
         *
         * @throws None
         */
        int n = this->no_subdivisions;
        return n*(k + n*(j + n*k));
    }

    int subdomain_id_to_linear_id(int i, int j, int k){
        /**
         * Converts a tuple of subdomain indices to a linear index.
         *
         * @param i the first subdomain index
         * @param j the second subdomain index
         * @param k the third subdomain index
         *
         * @return the linear index
         *
         * @throws None
         */
        int n = this->no_subdomains;
        return n*(k + n*(j + n*i)); 
    }


    std::array<int, 26> get_neighbouring_cells(int i, int j, int k){
        /**
         * Returns a length 26 integer array containing all neighbor cell ids.
         *
         * @param i the first cell index
         * @param j the second cell index
         * @param k the third cell index
         *
         * @return the length 26 integer array containing all neighbor cell ids
         *
         * @throws None
         */        
        std::array<int, 26> neighbouring_cells;
        // Take modulo of iu, and if id < 0 wrap back to no_subdivisions-1
        int iu = i+1 == this->no_subdivisions ? 0 : i+1;
        int id = i-1 == -1 ? this->no_subdivisions-1 : i-1;
        int ju = j+1 == this->no_subdivisions ? 0 : j+1;
        int jd = j-1 == -1 ? this->no_subdivisions-1 : j-1;
        int ku = k+1 == this->no_subdivisions ? 0 : k+1;
        int kd = k-1 == -1 ? this->no_subdivisions-1 : k-1;


        // Initialize all 26 combinations
        int idx = 0; 
        for (int di = -1; di <= 1; ++di) {
            for (int dj = -1; dj <= 1; ++dj) {
                for (int dk = -1; dk <= 1; ++dk) {
                    if (di == 0 && dj == 0 && dk == 0) continue; // Skip the current cell

                    // If an index is -1, set to no_subdivisions-1, if it is no_subdivisions, set to 0

                    int ni = (i + di + this->no_subdivisions) % this->no_subdivisions;
                    int nj = (j + dj + this->no_subdivisions) % this->no_subdivisions;
                    int nk = (k + dk + this->no_subdivisions) % this->no_subdivisions;

                    neighbouring_cells[idx++] = cell_id_to_linear_id(ni, nj, nk);
                }
            }
        }

        return neighbouring_cells;
    }
    

    private:
    // Variables, functions


}

#endif