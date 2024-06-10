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


class Communication{
    // Subdomains update arrays to be sent to other subdomains, kept in a communication object
    // To be used in the subdomain update function
    int no_subdomains; 

    // Communication steps are:
    // 1. Communicate ghost atoms; ids and positions (new included) of ghost atoms, we update all of them to avoid checking for ids


    // Data structures for forward and backward communication
    // Forward communication entails the migration of ghost atoms
    // Backward communication entails the reception of force calculation results
    // Keep an array of arrays of pointers to communication data stored in vectors of IDPOSITION structs
    std::vector<std::array<std::vector<IDPOSITION> , 26>> forward_communication_data;
    // TODO : std::vector<std::array<IDPOSITION *, 26>> backward_communication_data;

public:
    Communication(int no_subdomains) : no_subdomains(no_subdomains) {
        // Initialize the communication data structures
        forward_communication_data.resize(no_subdomains);
        backward_communication_data.resize(no_subdomains);

        for (int i = 0; i < no_subdomains; i++) {
            for (int j = 0; j < 26; j++) {
                forward_communication_data[i][j] = nullptr;
                backward_communication_data[i][j] = nullptr;
            }
        }
    }
    ~Communication() {
        // Destroy the communication data structures
        for (int i = 0; i < no_subdomains; i++) {
            for (int j = 0; j < 26; j++) {
                delete forward_communication_data[i][j];
                delete backward_communication_data[i][j];
            }
        }

        forward_communication_data.clear();
        backward_communication_data.clear();
    }

    // Functions
    void add_forward_communication_data(int subdomain_id, int cell_id, IDPOSITION *data) {
        forward_communication_data[subdomain_id][cell_id] = data;
    }
    void add_backward_communication_data(int subdomain_id, int cell_id, IDPOSITION *data) {
        backward_communication_data[subdomain_id][cell_id] = data;
    }

};


class ProcessorSubdomain{
    // TODO : particle information only in cell grid
    public:
        int processor_id;
        int subdomain_id;
        int no_particles; 
        int no_subdivisions; 
        int no_cells; 
        double xl; 
        double yl; 
        double zl; 
        double xu; 
        double yu; 
        double zu;
        double side_length;
        bool comm_forward; 
        Communication communication; // TODO add pointer to this shared object


        // Fixed parameters
        static const int no_contacts_per_particle = 6;  // For reservation of contact list size
        static const double dt;  // Time step 


        // Neighbouring subdomains
        std::array<int, 26> neighbouring_subdomains;
        std::array<std::array<int, 6>, 26> subdomain_cell_loop_limits;

        // Particle positions, velocities etc. 
        std::vector<std::array<double, 3>> positions;
        std::vector<std::array<double, 3>> velocities;

        // Neighbouring subdomains are ordered according to 
        // i-1, i+1, j-1, j+1, k-1, k+1, i-1j-1, i-1j+1, i-1k-1, i-1k+1, i+1j-1, i+1j+1, i+1k-1, i+1k+1


        // Cell lists 
        std::vector<int> cell_list_id;      // Cell id of each particle
        std::vector<std::array<int, 26>> neighbouring_cells;       // Each cell keeps the ids of its neighbours
        std::array<int , this->no_cells> cell_list_particle_no;    // Number of particles in each cell
        std::array<std::vector <int>, this->no_cells> cell_lists;  // Particles in each cell

        // Data structures:
        // Neighlists, cell grid, ghost atoms, ghost cells, ghost subdomains
        std::vector<std::tuple<int, int>> neighbour_list_info;     // Neighbour list of tuples (offset, number of neighbours)
        std::vector<int> neighbour_list;                            // Neighbour list of particle ids, indexed by neighbour_list_info

        ProcessorSubdomain(int processor_id, 
                            int subdomain_id, 
                            int no_particles, 
                            int no_cells, 
                            double xl, 
                            double yl, 
                            double zl, 
                            double xu, 
                            double yu, 
                            double zu, 
                            double side_length, 
                            std::array<int, 26> neighbouring_subdomains, 
                            Communication communication){
            // Initialize the processor subdomain
            this->processor_id = processor_id;
            this->subdomain_id = subdomain_id;
            this->no_particles = no_particles;
            this->no_subdivisions = no_subdivisions;
            this->no_cells = no_cells;
            this->xl = xl;
            this->yl = yl;
            this->zl = zl;
            this->xu = xu;
            this->yu = yu;
            this->zu = zu;
            this->side_length = side_length;
            this->neighbouring_subdomains = neighbouring_subdomains;
            this->cell_list_particle_no = std::array<int, this->no_cells>(0);
            this->communication = communication;

        }


    std::vector <int> get_cell_lists_boundary(int boundary_id){
        // Get the cell lists at the boundary according to the boundary number, with numbers ordered according to 
        std::vector <int> boundary_cell_lists;

        // Boundary sets start and end values of loop
        int istart = this->subdomain_cell_loop_limits[boundary_id][0];
        int iend = this->subdomain_cell_loop_limits[boundary_id][1];
        int jstart = this->subdomain_cell_loop_limits[boundary_id][2];
        int jend = this->subdomain_cell_loop_limits[boundary_id][3];
        int kstart = this->subdomain_cell_loop_limits[boundary_id][4];
        int kend = this->subdomain_cell_loop_limits[boundary_id][5];

        // Loop over the boundary sets
        for (int i = istart; i < iend; i++){
            for (int j = jstart; j < jend; j++){
                for (int k = kstart; k < kend; k++){
                    boundary_cell_lists.push_back(this->cell_lists[i][j][k]);
                }
            }
        }

        return boundary_cell_lists;

        
    }


    void set_boundary_loop_limits() {
        int n = this->no_subdivisions;
        this->subdomain_cell_loop_limits[0] = {0, 1, 0, n, 0, n};        // lower X
        this->subdomain_cell_loop_limits[1] = {n-1, n, 0, n, 0, n};      // upper X
        this->subdomain_cell_loop_limits[2] = {0, n, 0, 1, 0, n};        // lower Y
        this->subdomain_cell_loop_limits[3] = {0, n, n-1, n, 0, n};      // upper Y
        this->subdomain_cell_loop_limits[4] = {0, n, 0, n, 0, 1};        // lower Z
        this->subdomain_cell_loop_limits[5] = {0, n, 0, n, n-1, n};      // upper Z

        // Diagonal boundaries

        // Lower X, lower Y diagonals
        this->subdomain_cell_loop_limits[6] = {0, 1, 0, 1, 0, n};        // (0, 0, z)
        this->subdomain_cell_loop_limits[7] = {0, 1, n-1, n, 0, n};      // (0, n-1, z)
        this->subdomain_cell_loop_limits[8] = {0, 1, 0, n, 0, 1};        // (0, y, 0)
        this->subdomain_cell_loop_limits[9] = {0, 1, 0, n, n-1, n};      // (0, y, n-1)

        // Upper X, lower Y diagonals
        this->subdomain_cell_loop_limits[10] = {n-1, n, 0, 1, 0, n};     // (n-1, 0, z)
        this->subdomain_cell_loop_limits[11] = {n-1, n, n-1, n, 0, n};   // (n-1, n-1, z)
        this->subdomain_cell_loop_limits[12] = {n-1, n, 0, n, 0, 1};     // (n-1, y, 0)
        this->subdomain_cell_loop_limits[13] = {n-1, n, 0, n, n-1, n};   // (n-1, y, n-1)

        // Lower X, upper Y diagonals
        this->subdomain_cell_loop_limits[14] = {0, 1, n-1, n, 0, n};     // (0, n-1, z)
        this->subdomain_cell_loop_limits[15] = {0, 1, 0, 1, 0, n};       // (0, 0, z)
        this->subdomain_cell_loop_limits[16] = {0, 1, 0, n, 0, 1};       // (0, y, 0)
        this->subdomain_cell_loop_limits[17] = {0, 1, 0, n, n-1, n};     // (0, y, n-1)

        // Upper X, upper Y diagonals
        this->subdomain_cell_loop_limits[18] = {n-1, n, n-1, n, 0, n};   // (n-1, n-1, z)
        this->subdomain_cell_loop_limits[19] = {n-1, n, 0, 1, 0, n};     // (n-1, 0, z)
        this->subdomain_cell_loop_limits[20] = {n-1, n, 0, n, 0, 1};     // (n-1, y, 0)
        this->subdomain_cell_loop_limits[21] = {n-1, n, 0, n, n-1, n};   // (n-1, y, n-1)

        // Lower Z diagonals
        this->subdomain_cell_loop_limits[22] = {0, 1, 0, n, 0, 1};       // lower Z plane, lower X
        this->subdomain_cell_loop_limits[23] = {n-1, n, 0, n, 0, 1};     // lower Z plane, upper X

        // Upper Z diagonals
        this->subdomain_cell_loop_limits[24] = {0, 1, 0, n, n-1, n};     // upper Z plane, lower X
        this->subdomain_cell_loop_limits[25] = {n-1, n, 0, n, n-1, n};   // upper Z plane, upper X
}

    std::vector<IDPOSITION> gather_ghost_particles(int subdomain_id){
        std::vector<int> boundary_cell_lists = this->get_cell_lists_boundary(subdomain_id);
        no_cells = boundary_cell_lists.size();

        // Subdomain variables
        int no_ghost_particles = 0; 
        int idx = 0;

        // Loop over each cell adding them to the ghost list, but first calculate the number of particles to be communicated to preallocate the required array
        for (int j = 0; j < no_cells; j++){
            no_ghost_particles += this->cell_list_particle_no[boundary_cell_lists[j]];
        }
        // Preallocate the array of IDPOS pairs
        std::vector<IDPOSITION> ghost_particles(no_ghost_particles);
        for (int j = 0; j < no_cells; j++){
            for (int k = 0; k < this->cell_list_particle_no[boundary_cell_lists[j]]; k++){
                ghost_particles[idx].id = this->cell_list[boundary_cell_lists[j]][k];
                ghost_particles[idx].position = this->positions[this->cell_list[boundary_cell_lists[j]][k]];
                idx++;
            }
        }

        return ghost_particles;
        
    }

    void forward_communication(){
        // Forward communicate ghost atoms to surrounding subdomains
        // Loop over each neighbouring subdomain
        std::vector <int> boundary_cell_lists;
   
        for (int i = 0; i < 26; i++){
            std::vector<IDPOSITION> ghost_particles = this->gather_ghost_particles(i);
            communication[this->subdomain_id][this->neighbouring_subdomains[i]] = ghost_particles;
        }


    }

    void backward_communication(){
        // Backwards communicate force calculations on ghost atoms
    }




    private:

        void initial_init(){
            // Initialize particles
            this->init_particle_positions();
            this->init_particle_velocities();

            // Cell grid
            build_cell_lists();
            init_cell_grid();
            sort_by_cell_id();

        }






        // Particles initializations
        void init_particle_positions(){
            // Initialize the particles uniformly in the simulation box
            std::random_device rd;
            std::mt19937 gen(rd);
            std::uniform_real_distribution<> dis(0.0, side_length);

            for (int i = 0; i < this->no_particles; i++){
                this->positions[i][0] = xl + dis(gen);
                this->positions[i][1] = yl + dis(gen);
                this->positions[i][2] = zl + dis(gen);
            }
        }

        void init_particle_velocities(double avg_velocity){
            // TODO : Initialize the particles with velocity of avg_velocity in a random direction
            // We set these to zero for now

            for (int i = 0; i < this->no_particles; i++){
                this->velocities[i][0] = 0.0;
                this->velocities[i][1] = 0.0;
                this->velocities[i][2] = 0.0;

            }
        }



        // Cell grid initialization and management
        void build_cell_lists(){
            // Initialize the cell lists id for each particle, count number of particles in each cell and build cell lists
            this->cell_list_id = std::vector<int>(this->no_particles, 0);
            int lid = 0; 
            for (int i = 0; i < this->no_particles; i++){
                lid = this->pos_to_linear_cell_id(this->positions[i]);
                this->cell_list_particle_no[lid] += 1;
                this->cell_list_id[i] = lid;
                this->cell_lists[lid].push_back(i);
            }

        }
        // Neighbouring subdomains
        void init_neigbouring_subdomain(){
            // Initialize the neighbouring subdomain list of this subdomain
            
        }

        int pos_to_linear_cell_id(std::array<double, 3> pos){
            /**
             * Converts a position to a tuple of cell indices.
             *
             * @param pos the position to convert
             *
             * @return a tuple of three integers representing the cell indices (i, j, k)
             *
             * @throws None
             */
            // Proportions along each axis
            double px = pos[0]/this->side_length;
            double py = pos[0]/this->side_length;
            double pz = pos[0]/this->side_length;
                
            // These are mapped into the cells
            int i = int(px*this->no_subdivisions);
            int j = int(py*this->no_subdivisions);
            int k = int(pz*this->no_subdivisions);
            int n = this->no_cells;
            return n*(k + n*(j + n*i));
        }

        // Neighbor list initialization and management
        void init_neighbour_lists(){
            // Initialize the neighbour lists for each particle as a page containing the number of neighbours


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