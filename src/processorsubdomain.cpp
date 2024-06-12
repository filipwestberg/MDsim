#include <string>
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
#include <processorsubdomain.hpp>

// The configuration object keeping track of the current configuration
// of the program.

struct int3double
{
    int id;
    double position[3];
};

struct ij_pair
{
int i, j;
ij_pair(int i, int j) : i(i), j(j) {}
};

class ProcessorSubdomain{
    // TODO : particle information only in cell grid
    public:
        // Ids
        int processor_id;
        int subdomain_id;

        // Nos
        int no_particles; 
        int no_ghost_particles;
        int no_subdivisions; 
        int no_cells; 

        // Dimensions
        double xl; 
        double yl; 
        double zl; 
        double xu; 
        double yu; 
        double zu;

        double side_length;
        // Dimensions when including ghost cells
        double xlg; 
        double ylg; 
        double zlg; 
        double xug; 
        double yug; 
        double zug;

        double side_length_ghosts; 

        // Neighbour settings
        double neighbour_cutoff;
        double neighbour_cutoff_squared;

        bool comm_forward; 
        Communication communication; // TODO add pointer to this shared object

        std::function<double(std::array<double, 3>, std::array<double, 3>)> force_function;

        // Potential parameters
        std::vector<double> potential_parameters;

        // Fixed parameters
        static const int no_contacts_per_particle = 6;  // For reservation of contact list size
        static const double dt;  // Time step 
        static const double mass; // Mass of particle

        // Ghost particles: We increase the number of cells considered such that the cells composing the outermost shell are ghost cells in other domains, but not in the current domain.
        // When particles migrate into another domain, they migrate into the cells inside of the outermost shell. They are naturally incorporated into the subdomain.
        // Ghost particles can then efficiently be included during the building of the neighbour lists.
    
        // Global particle ids
        std::vector<int> global_particle_ids;

        // Neighbouring subdomains
        std::array<int, 26> neighbouring_subdomains;
        std::array<std::array<int, 6>, 26> subdomain_cell_loop_limits;

        // Particle positions, velocities etc. passed from the overarching configuration object
        std::vector<std::array<double, 3>> positions;
        std::vector<std::array<double, 3>> velocities;

        // Ghost particles
        std::vector<std::array<double, 3>> ghost_positions;
        
        // Cell lists 
        std::vector<int> cell_list_id;      // Cell id of each particle
        std::vector<std::array<int, 26>> neighbouring_cells;       // Each cell keeps the ids of its neighbours
        std::array<int , this->no_cells> cell_list_particle_no;    // Number of particles in each cell
        std::array<std::vector <int>, this->no_cells> cell_lists;  // Particles in each cell

        // Neighlists, cell grid, ghost atoms, ghost cells, ghost subdomains
        std::vector<std::tuple<int, int>> neighbour_list_info;     // Neighbour list of tuples (offset, number of neighbours)
        std::vector<ij_pair> neighbour_list;                       // Neighbour list of particle id pairs indicating neighborship

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

            // Set class variables
            this->processor_id = processor_id;
            this->subdomain_id = subdomain_id;
            this->no_particles = no_particles;
            this->no_ghost_particles = 0;
            this->no_subdivisions = no_subdivisions + 2;
            this->no_cells = (this-> no_subdivisions)^3;

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

            // Initialize structures, TODO : move to constructor, look over
            this->cell_lists = std::array<std::vector <int>;
            this->neighbouring_cells = std::array<std::array<int, 26>;
            this->no_cells>(std::vector <int> (0));
            this->neighbour_list = std::vector<ij_pair>(0);
            this->neighbour_list_info = std::vector<std::tuple<int, int>>(0);
            this->cell_list_id = std::vector<int>(0);
            this->subdomain_cell_loop_limits = std::array<std::array<int, 6>, 26>(std::array<int, 6>(0));   
            this->cell_lists = std::array<std::vector<int>, this->no_cells>(std::vector<int>(0));
            this->neighbouring_cells = std::array<std::array<int, 26>, this->no_cells>(std::array<int, 26>(0));


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
            this->communication.set_forward_communication_data(this->subdomain_id, this->neighbouring_subdomains[i],  ghost_particles)
        }


    }

    void backward_communication(){
        // Backwards communicate force calculations on ghost atoms
    }




    private:

        void initial_init(){
            // Cell grid
            init_neighbouring_cell_lists();
            build_cell_lists();
            build_neighbouring_cell_lists();
            init_cell_grid();
            sort_by_cell_id();

            // Neighbour list
            init_neighbour_list();
            build_neighbour_lists();

        }




        // ---------------------------------- Cell grid initialization and management --------------------------------

        void init_neighbouring_cell_list(){
            // Initialize the empty neighbour list
            this->neighbour_list.clear();
        }
        void init_cell_lists(){
            // Initialize the empty neighbour lists for each cell
            this->cell_lists = std::vector<std::vector<int>>(this->no_cells, std::vector<int>());
            this->cell_list_id = std::vector<int>(this->no_particles, 0);
            this->cell_list_particle_no = std::vector<int>(this->no_cells, 0);
        }
        void build_neighbouring_cell_lists(){
            // Initialize the neighbouring cell lists of increasing i, j, k (templating to reduce memory usage)
            std::tuple<int, int, int> current_cell_id;
            int lid; 

            for (int i = 0; i < this->no_cells; i++){
                current_cell_id = this->linear_cell_id_to_pos(i);
                for (int j = 0; j < 1; j++){
                    for (int k = 0; k < 1; k++){
                        for (int l = 0; l < 1; l++){
                            if (j == 0 && k == 0 && l == 0){
                                continue;
                            } else{
                                lid = this->pos_to_linear_cell_id({std::get<0>(current_cell_id) + j, std::get<1>(current_cell_id) + k, std::get<2>(current_cell_id) + l});
                                this->neighbouring_cell_lists[i].push_back(this->pos_to_linear_cell_id({std::get<0>(cell_id) + j, std::get<1>(cell_id) + k, std::get<2>(cell_id) + l}));

                            }
                        }
                    }
                    
                }
            }

        }

        void build_cell_lists(){
            // Initialize the cell lists id for each particle, count number of particles in each cell and build cell lists
            // Initialize an empty int vector
            
            int lid = 0; 
            for (int i = 0; i < this->no_particles; i++){
                lid = this->pos_to_linear_cell_id(this->positions[i]);
                this->cell_list_particle_no[lid] += 1;
                this->cell_list_id[i] = lid;
                this->cell_lists[lid].push_back(i);
            }

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
            if (pos[0] < this->xgl || pos[1] < this->ylg || pos[2] < this->zlg){
                std::cerr << "Position out of bounds" << std::endl;
                exit(1);
            }

            // Proportions along each axis in [xl-sidelength/no_subdivisions, xu+sidelength/no_subdivisions]
            double px = (pos[0]-this->xgl)/this->side_length_ghosts;
            double py = (pos[1]-this->ylg)/this->side_length_ghosts;
            double pz = (pos[2]-this->zlg)/this->side_length_ghosts;
                
            // These are mapped into the cells (remember that the ghost cells are included as the outermost cells)
            int i = int(px*this->no_subdivisions);
            int j = int(py*this->no_subdivisions);
            int k = int(pz*this->no_subdivisions);
            int n = this->no_cells;

            return n*(k + n*(j + n*i));
        }

        // ---------------------------------- Neighbor list initialization and management --------------------------------
        void init_neighbour_list(){
            // Initialize the empty neighbour list
            this->neighbour_list = std::vector<ij_pair>(this->no_particles, 0);
        }
        
        void check_neighbours(const std::vector<int>& particle_ids1, const std::vector<int>& particle_ids2){
            // Given the ids in the cell lists of two particles, check if there are any neighbours and if so add them to the neighbour list
            std::array<double, 3> pos1;
            std::array<double, 3> pos2;
            double dist_squared;
            // Iterate over each particle in the array
            for (int i = 0; i < particle_ids1.size(); i++){
                pos1 = this->positions[particle_ids1[i]];
                // Iterate over each particle in the other array
                for (int j = 0; j < particle_ids2.size(); j++){
                    pos2 = this->positions[particle_ids2[j]];
                    dist_squared = (pos1[0] - pos2[0])*(pos1[0] - pos2[0]) + 
                                   (pos1[1] - pos2[1])*(pos1[1] - pos2[1]) + 
                                   (pos1[2] - pos2[2])*(pos1[2] - pos2[2]);

                    if (dist_squared < this->neighbour_cutoff_squared){
                        // Add the ijpair to the neighbour list
                        this->neighbour_lists[i].push_back(ij_pair({i, j}));
                    }
                }
            }
            
        }

        void build_neighbour_lists(){
            // Build the neighbour lists for each particle, based on in and out of subdomain particles
            // Iterate over each cell, and then over each particle in the cell
            for (int i = 0; i < this->no_cells; i++){
                // In cell particle neighbour
                check_neighbours(this->cell_lists[i], this->cell_lists[i]);
                // Out of cell neighbours, iterate over the neighbouring cells that are of increasing i, j, k
                for (int j : neighbouring_cell_lists[i]){
                    check_neighbours(this->cell_lists[neighbouring_cell_lists[i][j]], this->cell_lists[i]);
                }
            }
        }


        //----------------------------------- Force calculations -----------------------------------
        void force_sweep(){
            // Calculate forces on all particles in the subdomain
            // Iterate over the neighbour list
            // Zero the forces initially
            for (int i = 0; i < this->no_particles; i++){
                this->forces[i][0] = 0.0;
                this->forces[i][1] = 0.0;
                this->forces[i][2] = 0.0;
            }

            int p1; 
            int p2; 
            std::array<double, 3> pos1;
            std::array<double, 3> pos2;
            std::array<double, 3> force;
            // Sweep over the subdomain neighbour list (including ghost particles)
            for (int i = 0; i < this->neighbour_lists.size(); i++){
                p1 = this->neighbour_lists[i][0].i;
                p2 = this->neighbour_lists[i][0].j;
                pos1 = this->positions[p1];
                pos2 = this->positions[p2];
                force = this->force_function(pos1, pos2);

                // Add the force and its opposite to the other particle
                this->forces[p1][0] += force[0];
                this->forces[p1][1] += force[1];
                this->forces[p1][2] += force[2];

                this->forces[p2][0] -= force[0];
                this->forces[p2][1] -= force[1];
                this->forces[p2][2] -= force[2];
                
            }
        }


        double force_function(std::array<double, 3> pos1, std::array<double, 3> pos2){
            if (this->force_function){
                return this->force_function(pos1, pos2);
            } else {
                std::cerr << "Force function not set" << std::endl;
                return 0.0;
            }
            
        }

        void set_force_parameters(std::vector<double> force_parameters, const std::string& force_function_name){
            // Set the parameters of the force function
        }

        void set_force_function(const std::string& force_function_name) {
                // Define the map of force functions
                const std::unordered_map<std::string, ForceFunction> force_function_map{
                    {"harmonic", harmonic},
                    {"LennardJones", LennardJones},
                    {"Hooke", Hooke}
                };

                // Find the force function by name
                const auto it = force_function_map.find(force_function_name);
                if (it != force_function_map.end()) {
                    this->force_function = it->second;
                } else {
                    std::cerr << "Force function not found: " << force_function_name << std::endl;
                    this->force_function = nullptr;  // Or some default behavior
                }
            }
            
        void set_force_arguments(std::vector<double> force_arguments){
            // Set and parse the arguments of the force function, TODO


        }
        

        // TODO: fix argument parsing and force functions
        double harmonic(std::array<double, 3> pos1, std::array<double, 3> pos2){
            return 0.0;
        }

        double LennardJones(std::array<double, 3> pos1, std::array<double, 3> pos2){
            return 0.0;
        }

        double Hooke(std::array<double, 3> pos1, std::array<double, 3> pos2){
            return 0.0;
        }

        // ------------------------------ Time integration -------------------------------
        void update_positions(double dt){
            // Verlet integrate the positions
            for (int i = 0; i < this->no_particles; i++){
                this->positions[i][0] += this->velocities[i][0] * dt + 0.5 * this->forces[i][0]/this->mass * dt * dt;
                this->positions[i][1] += this->velocities[i][1] * dt + 0.5 * this->forces[i][1]/this->mass * dt * dt;
                this->positions[i][1] += this->velocities[i][2] * dt + 0.5 * this->forces[i][2]/this->mass * dt * dt;
            }
        }
        // ------------------------------ Communication -------------------------------
        void unpack_ghost_particles(){
            // Unpack a subdomain's ghost particles into the cell lists
        }

        // ------------------------------ Auxillary functions -------------------------------
        void particle_no_tracking() {
            // Track the number of particles and ghost particles in this subdomain
            int no_ghost_particles = 0;
            int n = this->no_subdivisions;

            // Iterate over the ghost cells and add up all the ghost particles
            // Ghost cells are at the edges of the simulation box

            // Iterate over the cells on the i-plane (along z-axis, x and y vary)
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    no_ghost_particles += this->cell_list_particle_no[this->pos_to_linear_cell_id({i, j, 0})];
                    no_ghost_particles += this->cell_list_particle_no[this->pos_to_linear_cell_id({i, j, n-1})];
                }
            }

            // Iterate over the cells on the j-plane (along y-axis, x and z vary)
            for (int i = 0; i < n; i++) {
                for (int k = 1; k < n-1; k++) { // Skip the corners already counted in the i-plane loop
                    no_ghost_particles += this->cell_list_particle_no[this->pos_to_linear_cell_id({i, 0, k})];
                    no_ghost_particles += this->cell_list_particle_no[this->pos_to_linear_cell_id({i, n-1, k})];
                }
            }

            // Iterate over the cells on the k-plane (along x-axis, y and z vary)
            for (int j = 1; j < n-1; j++) { // Skip the corners already counted in the i-plane loop
                for (int k = 1; k < n-1; k++) { // Skip the edges already counted in the j-plane loop
                    no_ghost_particles += this->cell_list_particle_no[this->pos_to_linear_cell_id({0, j, k})];
                    no_ghost_particles += this->cell_list_particle_no[this->pos_to_linear_cell_id({n-1, j, k})];
                }
            }


            // Set these
            this->no_ghost_particles = no_ghost_particles;
            this->no_particles = this->positions.size() - this->no_ghost_particles;


        }


}

