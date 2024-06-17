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
#include <unordered_map>
#include <unordered_set>
#include <omp.h>
#include <map>
#include <processorsubdomain.hpp>

// The configuration object keeping track of the current configuration
// of the program.
// TODO : 
//          Move force functions and binding etc to configuration and pass a lambda to this class
//          Periodicity of the system, movement and neighborship
//          Coherency of all structures (id, ghost_comm_array etc).
//          Ghost cutoff as variable
//          Identify particles shifting subdomain
//          Incorporate particles shifting subdomain
//          Exception handling
//          Delete particle : shift the particles -> delete last element in the array to avoid having to change all local ids

struct int3double    
{
    int id;
    double position[3];
    int3double(int id, double position[3]) : id(id), position(position) {}
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
        int no_local_particles; 
        int no_ghost_particles;
        int no_ghost_particles_others; 
        int no_subdivisions; 
        int no_subdivisions_ghosts; 
        int no_cells; 
        int no_cells_ghosts;

        // Dimensions
        double bounds[6];
        double gh_bounds[6];

        double cell_length;
        double side_length;

        // Midpoint
        double xmid;
        double ymid;
        double zmid;

        // Neighbour settings
        double neighbour_cutoff;
        double neighbour_cutoff_squared;

        // Misc.
        int iteration = 0; // Current iteration
        double ghost_cutoff; // Cutoff for ghost particles
        int  neighbor_list_build_per = 5;

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
    
        // Global particle ids for in domain particles and ghost particles
        std::map<int, int> global_local_particle_ids_map;
        std::map<int, int> local_global_particle_ids_map;
        std::map<int, int> global_local_ghost_particle_ids_map;
        std::map<int, int> local_global_ghost_particle_ids_map;

        std::vector<int> global_particle_ids_vector;        // Initial global particle ids for in domain particles

        // Neighbouring subdomains
        std::array<int, 6> neighbouring_subdomains;

        // Particle positions, velocities etc. passed from the overarching configuration object
        std::vector<std::array<double, 3>> positions;
        std::vector<std::array<double, 3>> velocities;

        // Cell lists 
        std::vector<int> cell_list_id;                             // Cell id of each particle
        std::vector<std::array<int, 26>> neighbouring_cells;       // Each cell keeps the ids of its neighbours
        std::array<int , this->no_cells> cell_list_particle_no;    // Number of particles in each cell
        std::array<std::vector <int>, this->no_cells> cell_lists;  // Particles in each cell
        std::array<std::vector<std::array<int, 3>>, this->no_cells> cell_list_pos; // Particle positions in each cell

        // Neighlists, cell grid, ghost atoms, ghost cells, ghost subdomains
        std::vector<ij_pair> neighbour_list;                       // Neighbour list of particle id pairs indicating neighborship
        std::vector<ij_pair> ghost_neighbour_list;                // Neighbour list of ghost particle id pairs indicating neighborship

        // Movement array
        std::vector<double> accumulated_movement;

        // Communication arrays 
        std::array<std::vector<int3double> , 6> forward_communication_data;
        std::array<std::vector<int3double> , 6> backward_communication_data;
        std::vector<bool> ghost_comm_array;

        // Flags
        bool initialized = false;
        bool dump_data = false; 
        bool rebuild_neighbours = false;

        // Removal queue
        std::queue<int> removal_queue;

        // Reordering queue (local particle ids that are to be moved either into or out of the ghost particle subarray)
        std::vector<int> reordering_queue;

        // Ghost addition queue and global id set
        std::queue<int> ghost_addition_queue_id;
        std::queue<double[3]> ghost_addition_queue_pos;
        std::map<int, bool> ghost_removal_map;

        // Position array reordering
        std::vector<int> ghost_reordering_ids;
        std::vector<int> in_domain_reordering_ids;


        ProcessorSubdomain(int processor_id, 
                            int subdomain_id, 
                            int no_particles, 
                            int no_cells, 
                            double bounds[6], 
                            double side_length, 
                            std::array<int, 26> neighbouring_subdomains,
                             
                            Communication communication) : this->initializations(){

            // Set class variables
            this->processor_id = processor_id;
            this->subdomain_id = subdomain_id;
            this->no_particles = no_particles;
            this->no_ghost_particles = 0;

            this->no_subdivisions = no_subdivisions;
            this->no_subdivisions_ghosts = no_subdivisions + 2;
            this->no_cells = std::pow(this->no_subdivisions, 3);
            this->no_cells_ghosts = std::pow(this->no_subdivisions_ghosts, 3);

            this->bounds = bounds;

            this->side_length = side_length;
            this->neighbouring_subdomains = neighbouring_subdomains;
            this->communication = communication;
        }

    private:
        // ------------------------------- Iteration ------------------------------------------------
        // Over the previous timestep already known ghost particles have moved, and have either been introduced or removed from the set of ghost particles.
        // In domain particles have moved and have either been introduced or removed from the set of ghost particles in other domains.
        // In domain particles may have moved into another domain.
        // 
        // Iterate the subdomain one step by : 
        // 1. Unpacking ghost particles:
        //                               1.a) add ghost particles onto the ghost positions array
        //                               1.b) move ghost particles
        //                               1.b) remove ghost particles from the ghost particle array that have not been communicated this timestep
        //                               1.d) particles that have moved across the boundary are to be added into the local position array 
        // 2. Check the integrity of the particle positions array
        // 3. Rebuild the neighbouring lists if necessary
        // 4. Force calculation
        // 5. Integrate
        // 6. Communicate forward ghost particles (this includes genuine ghost particles and particles that have moved across the boundary)



        void iterate(){
            // Unpack ghost particles
            this->unpack_ghost_particles();

            // Check the integrity of the particle positions array and reorder if necessary (former ghost particles are here made into in-domain particles)
            this->order_positions();

            // Rebuild neighbour lists on nth iteration
            if (this->iteration % this->neighbor_list_build_period == 0) {this->build_neighbour_lists();}
            
            // Force calculation
            this->force_sweep();

            // Integrate
            this->update_positions();

            // Communicate forward ghost particles
            this->communicate_ghost_particles();

        }
        
        // ---------------------------------- Initialization --------------------------------
        {

        void init_particle_mapping(){
            // Initialize the global -> local and local -> global particle mappings from the global particle id vector given at construction
            for (int i = 0; i < this->no_particles; i++){
                this->global_local_particle_ids_map[this->global_particle_ids_vector[i]] = i;
                this->local_global_particle_ids_map[i] = this->global_particle_ids_vector[i];
            }
        }
        }
    
        // ---------------------------------- Cell grid initialization and management --------------------------------
        {
        // Initialize the :
        //                     - empty neighbour list
        //                     - empty cell lists
        //                     - empty positions array
        //                     - empty communication arrays
        //                     - empty neighbouring cell lists

        void initializations(){
            // Initialize the empty forward and backward communication arrays
            for (int i = 0; i < 6; i++){
                this->forward_communication_data[i] = std::vector<int3double>();
                this->backward_communication_data[i] = std::vector<int3double>();
            }

            // Initialize the empty ghost communication array
            this->ghost_comm_array = std::vector<bool>;

            // Initialize the empty neighbour list
            this->neighbour_list = std::vector<ij_pair>;
            this->ghost_neighbour_list = std::vector<ij_pair>;

            // Initialize the empty cell lists, the cell positions, neighbouring cells 
            this->cell_lists = std::vector<std::vector<int>>(this->no_cells_ghosts, std::vector<int>());
            this->cell_list_particle_no = std::vector<int>;        

            // Initialize the empty positions and velocity arrays of size no_particles
            this->positions = std::vector<std::array<double, 3>>(this->no_particles, {0.0, 0.0, 0.0});
            this->velocities = std::vector<std::array<double, 3>>(this->no_particles, {0.0, 0.0, 0.0});

            // Initialize the local - global particle mappings
            this->local_global_particle_ids_map = std::map<int, int>();
            this->global_local_particle_ids_map = std::map<int, int>();

            this->global_local_ghost_particle_ids_map = std::map<int, int>();
            this->local_global_ghost_particle_ids_map = std::map<int, int>();

            // Initialize the neighbouring subdomains
            this->neighbouring_subdomains = std::array<int, 6>(0);

            // Initialize the removal queue
            this->removal_queue = std::queue<int>;

        }

        void parameter_setup(){
            // Setup midpoints
            this->midpoint = std::array<double, 3>({this->bounds[0] + this->side_length/2.0, this->bounds[2] + this->side_length/2.0, this->bounds[4] + this->side_length/2.0});

        }


        void build_cell_lists(){
            // Initialize the cell lists id for each particle, count number of particles in each cell and build cell lists
            // Include in domain particles
            int lid = 0; 
            for (int i = 0; i < this->no_particles; i++){
                lid = this->pos_to_linear_cell_id(this->positions[i]);
                this->cell_list_particle_no[lid] += 1;
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
                std::cerr << "Position out of bounds of subdomain, including ghost cells" << std::endl;
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

            return this->no_cells*(k + this->no_cells*(j + this->no_cells*i));
        }
        }
        // ---------------------------------- Neighbor list initialization and management --------------------------------
        {
        void check_add_neighbours(const std::vector<int>& cell_ids1, 
                              const std::vector<int>& cell_ids2, 
                              const std::vector<std::array<int, 3>>& positions1,
                              const std::vector<std::array<int, 3>>& positions2){
            // Given the ids in the cell lists of two particles, check if there are any neighbours and if so add them to the neighbour list
            std::array<double, 3> pos1;
            std::array<double, 3> pos2;
            double dist_squared;
            // Iterate over each particle in the array
            for (int i = 0; i < cell_ids1.size(); i++){
                pos1 = this->positions1[cell_ids1[i]];
                // Iterate over each particle in the other array
                for (int j = 0; j < cell_ids2.size(); j++){
                    pos2 = this->positions2[cell_ids2[j]];
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
            int j;
            for (int i = 0; i < this->no_cells; i++){
                // In cell particle neighbour
                check_add_neighbours(this->cell_lists[i], this->cell_lists[i]);

                // Out of cell neighbours, iterate over the neighbouring cells
                // Make use of a half template
                j = ijk_to_lid(i+1, j, k);
                check_add_neighbours(this->cell_lists[neighbouring_cell_lists[i][j]], this->cell_lists[i], this->positions, this->positions);

                j = ijk_to_lid(i, j+1, k);
                check_add_neighbours(this->cell_lists[neighbouring_cell_lists[i][j]], this->cell_lists[i], this->positions, this->positions);

                j = ijk_to_lid(i, j, k+1);
                check_add_neighbours(this->cell_lists[neighbouring_cell_lists[i][j]], this->cell_lists[i], this->positions, this->positions);

                j = ijk_to_lid(i+1, j+1, k);
                check_add_neighbours(this->cell_lists[neighbouring_cell_lists[i][j]], this->cell_lists[i], this->positions, this->positions);

                j = ijk_to_lid(i+1, j, k+1);
                check_add_neighbours(this->cell_lists[neighbouring_cell_lists[i][j]], this->cell_lists[i], this->positions, this->positions);

                j = ijk_to_lid(i, j+1, k+1);
                check_add_neighbours(this->cell_lists[neighbouring_cell_lists[i][j]], this->cell_lists[i], this->positions, this->positions);

                j = ijk_to_lid(i+1, j+1, k+1);
                check_add_neighbours(this->cell_lists[neighbouring_cell_lists[i][j]], this->cell_lists[i], this->positions, this->positions);
                
            }
        }

        int ijk_to_lid(int i, int j, int k){
            return this->no_cells*(k + this->no_cells*(j + this->no_cells*i));
        }

        }
        // ----------------------------------- Ghost particles -----------------------------------
        {




        }

        //----------------------------------- Force calculations -----------------------------------
        {
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
                if (p1<this->no_local_particles){
                    this->forces[p1][0] += force[0];
                    this->forces[p1][1] += force[1];
                    this->forces[p1][2] += force[2];
                }

                if (p2<this->no_local_particles){
                    this->forces[p2][0] -= force[0];
                    this->forces[p2][1] -= force[1];
                    this->forces[p2][2] -= force[2];
                }
                
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
        }
        // ------------------------------ Time integration -------------------------------
        {
        void update_positions(double dt){
            // Verlet integrate the positions
            for (int i = 0; i < this->no_local_particles; i++){
                this->positions[i][0] += this->velocities[i][0] * dt + 0.5 * this->forces[i][0]/this->mass * dt * dt;
                this->positions[i][1] += this->velocities[i][1] * dt + 0.5 * this->forces[i][1]/this->mass * dt * dt;
                this->positions[i][1] += this->velocities[i][2] * dt + 0.5 * this->forces[i][2]/this->mass * dt * dt;
            }
        }
        }
        // Particle removal, addition, editing

        void add_ghost_particles(){
            // Add the ghost particles to the position array, remember that they *can* be inside the domain. We fix this at a later point in the iteration.
            while(!this->ghost_addition_queue.empty()){
                // Add the ghost particles
                const auto global_id = this->ghost_addition_queue.front().first;
                const auto& position = this->ghost_addition_queue.front().second;
                this->add_particle(global_id, position);
                this->ghost_addition_queue.pop_front();
            }
        }

        void remove_ghost_particles(){
            // Remove the ghost particles from the position array based on the removal map
            for (const auto & [global_id, state] : this->ghost_removal_map) {
                if (state) {
                    _remove_ghost_particle_(global_id);
                }
            }
        }
        void remove_particle(int global_id, bool ghost = false){
            // Remove the given particle from the subdomain
            if (!ghost) {
                _remove_domain_particle_(global_id);
            } else {
                _remove_ghost_particle_(global_id);
            }

        }
        void add_particle(int global_id, std::array<double, 3> position, bool ghost = false){
            // Add the given particle to the subdomain
            if (!ghost) {
                _add_domain_particle_(global_id, position);
            } else {
                _add_ghost_particle_(global_id, position);
            }

        }
        void _add_ghost_particle_(int global_id, std::array<double, 3> position){
            // Add the given ghost particle to the subdomain
            this->positions.push_back(position);
            // New local id 
            int local_id = this->positions.size()-1;
            this->global_local_ghost_particle_ids_map[global_id] = local_id;
            this->local_global_ghost_particle_ids_map[local_id] = global_id;
            // Add it to the removal map
            this->ghost_removal_map[local_id] = false; 
        }
        void _add_domain_particle_(int global_id, std::array<double, 3> position){
            // Add the given particle to the subdomain
            this->positions.push_back(position);
            // New local id 
            int local_id = this->positions.size()-1;
            this->global_local_particle_ids_map[global_id] = local_id;
            this->local_global_particle_ids_map[local_id] = global_id;

        }

        void remove_particle(int global_id){
            // Remove the given particle from the subdomain
            int local_id = this->global_local_particle_ids_map[global_id];
            this->positions.erase(this->positions.begin() + local_id);
            this->global_local_particle_ids_map.erase(global_id);
            this->local_global_particle_ids_map.erase(local_id);
            
        }


        void order_positions(){
            // Fix the order of the positions in the position array
            verify_positions();
            reorder_positions();
        }

        void verify_positions(){
            // Check the order of the particles in the position array
            bool in_domain;
            bool is_ghost;
            // Verify size of the position array
            veriify_size();

            // First supposed in-domain particles
            for (int i = 0; i < this->no_local_particles; i++){
                in_domain = this->inside_domain(this->positions[i]);
                if (!in_domain) {
                    is_ghost = is_ghost_particle(this->positions[i]);
                    if (is_ghost) {
                        // This is a ghost particle
                        this->in_domain_reordering_ids.push_back(i);
                    } else {
                        // This particle should be deleted; add it to the removal queue
                        this->removal_ids.push_back(i);
                    }
                    
                } else {
                    continue;
                }
            }
            // Then supposed ghost particles
            for (int i = this->no_local_particles; i < this->no_particles; i++){
                is_ghost = is_ghost_particle(this->positions[i]);
                if (!is_ghost) {
                    // Ghost particle not as it should  be
                    in_domain = this->inside_domain(this->positions[i]);
                    if (in_domain) {
                        // A local particle! 
                        this->ghost_reordering_ids.push_back(i);
                    } else {
                        // This particle should be removed
                        this->removal_ids.push_back(i);
                    }
                
                } else {
                continue;
                }
            }
        }


        void reorder_positions(){
            // Based on the local and ghost reordering queues reorder the position array
            // Edge cases: TODO
            int no_local = this->in_domain_reordering_ids.size();
            int no_ghost = this->ghost_reordering_ids.size();
            if (no_ghost == 0 && no_local == 0) {
                // No reordering
                return;
            } else if (no_local > 0 && no_ghost == 0) {
                // Only local reordering, exchange the position of the last local particle with this one, subtract one from the number of local particles and continue
                exchange_position(global_id_a = this->in_domain_reordering_ids.front(), global_id_b = this->ghost_reordering_ids.front(), local_a = true, local_b = false);
                this->in_domain_reordering_ids.pop_front();
                this->no_local_particles -= 1;
            } else if (no_ghost > 0 && no_local == 0) {
                // Only ghost reordering, exchange the position of the first ghost particle with this one, increase the local particle counter and continue
                exchange_position(this->ghost_reordering_ids.front(), this->no_local_particles);
                exchange_position(global_id_a = this->ghost_reordering_ids.front(), global_id_b = this->no_local_particles, local_a = false, local_b = true);
                this->ghost_reordering_ids.pop_front();
                this->no_local_particles += 1;
            } else {
                // Do one one exchange between the queues and recursively call this function, TODO : make this more efficient through explicit code
                exchange_position(this->in_domain_reordering_ids.front(), this->ghost_reordering_ids.front());
                exchange_position(global_id_a = this->in_domain_reordering_ids.front(), global_id_b = this->ghost_reordering_ids.front(), local_a = true, local_b = false);
                this->in_domain_reordering_ids.pop_front();
                this->ghost_reordering_ids.pop_front();
                reorder_positions();
            }
            return;
        }

        void exchange_position(int global_id_a, 
                                int global_id_b, 
                                bool local_a = false, 
                                bool local_b = false){
                                            /**
             * Exchange the positions of two particles.
             *
             * @param particle_a the id of the first particle
             * @param particle_b the id of the second particle
             * @param local_a flag indicating if a is local, if not its assumed to be a ghost
             * @param local_b flag indicating if b is local, if not its assumed to be a ghost
             */

            // Exchange the positions of two particles
            std::array<double, 3> position_a = this->positions[global_id_a];
            std::array<double, 3> position_b = this->positions[global_id_b];
            this->positions[global_id_a] = position_b;
            this->positions[global_id_b] = position_a;

            // We now need to update the maps based on the state change of the particle
            local_id_a = (local_a ? this->global_local_particle_ids_map : this->global_local_ghost_particle_ids_map)[global_id_a];
            local_id_b = (local_b ? this->global_local_particle_ids_map : this->global_local_ghost_particle_ids_map)[global_id_b];

            (local_a ? this->local_global_particle_ids_map : this->local_global_ghost_particle_ids_map)[local_id_a] = global_id_b;
            (local_b ? this->local_global_particle_ids_map : this->local_global_ghost_particle_ids_map)[local_id_a] = global_id_a;

            (local_a ? this->global_local_particle_ids_map : this->global_local_ghost_particle_ids_map)[global_id_a] = local_id_a;
            (local_b ? this->global_local_particle_ids_map : this->global_local_ghost_particle_ids_map)[global_id_b] = local_id_b;

            
            }

        void verify_size(){
            // Check that no_ghosts + no_locals = no_particles
            bool test1 = this->no_ghosts + this->no_locals == this->no_particles;
            bool test2 = this->no_ghosts + this->no_locals == this->positions.size();
            if (!test1 || !test2) {
                // Something has gone wrong
                throw std::runtime_error("ProcessorSubdomain: Size verification failed");
            }
        }
        // ------------------------------ Map Interaction -------------------------------
        {
            exchange_map_ids(int global_id_a, 
                             int global_id_b, 
                             auto& map_a_global_local, 
                             auto& map_b_global_local, 
                             auto& map_a_local_global, 
                             auto& map_b_local_global) {
            // TODO : Create a two-way associative map (local - global) to (global - local)
            // Exchange local ids in the maps
            local_id_a = map_a_global_local[global_id_a];
            local_id_b = map_b_global_local[global_id_b];
            // Exchange local ids in the global map
            map_a_local_global[local_id_a] = global_id_b;
            map_b_local_global[local_id_b] = global_id_a;
            // Exchange global ids in the maps
            map_a_global_local[global_id_b] = local_id_a;
            map_b_global_local[global_id_a] = local_id_b;
        }
        }

        // ------------------------------ Communication -------------------------------
        {

        void prepare_ghost_particles() {
            // Set the map state of the ghosts, such that we can later remove the ghost particles that have migrated out of the ghost zones
            set_map_state(false, this->ghost_removal_map);
            // Unpack the ghost particles, editing particle positions
            unpack_ghost_particles();
            // Add the ghost particles in the ghost buffer
            add_ghost_particles();
            // Remove the ghost particles that have migrated out of the ghost zones
            remove_ghost_particles();
        }

  
        void set_map_state(bool map_state, std::map<auto, bool> map) {
            // Iterate over all elements in the map and set their value to the given state
            for (auto& [key, value] : map) {
                // Set the value of the current element to the given state
                value = map_state;
            }
        }

        // Ghost lambda functions in x, y, z directions
        auto chghost_dir = [] (double r, double l, double u) {return (r > l && r < u); };

        void unpack_ghost_particles() {
            // Unpack a subdomain's ghost particles: editing particle positions, adding ghost particles and removing ghost particles
            // Messages are vectors of IDPOSITION structs
            // Iterate over neighbouring subdomains

            // Iterate over the neighboring subdomains
            for (int j : this->neighbouring_subdomains) {
                // Get the communication data for the current neighboring subdomain
                // Check if there are any ghost particles to unpack

                if (this->communication.get_forward_communication_data_size(this->subdomain_id, j) == 0) continue;

                const auto& comm_data = this->communication.get_forward_communication_data(this->subdomain_id, j);

                
                // Iterate over the ghost particles in the communication data
                for (const auto& particle : comm_data) {
                    int global_id = particle.first;
                    const auto& position = particle.second;

                    // Check if the global ID is in the global-to-local ghost particle ID map
                    auto containsID = this->global_local_ghost_particle_ids_map.contains(global_id);
                    if (containsID) {
                        // Particle exists, update its position
                        int local_id = it->second;
                        this->positions[local_id] = position;
                        this->ghost_removal_map[local_id] = false; 
                    } else {
                        // New particle, buffer it into the ghost addition queue
                        this->ghost_addition_queue_ids.push_back(global_id);
                        this->ghost_addition_queue_positions.push_back(position);
                    }
                }
                
            }
        }

        
        void pack_ghost_particles_X(){
            // Pack this subdomain's ghost particles into messages to its neighbouring subdomains in the X direction
            // Accumulate the ghost particles in the communication arrays
            int3double idpos;
            bool lx, ux;
            for (int i = 0; i < this->no_local_particles; i++) {
                // Check for ghost particles in each direction and pack them accordingly. 
                // Particles that have moved into another subdomain are communicated here   

                idpos = int3double(i, this->positions[i]);

                // Check if the particle is across x-boundaries
                lx = position[i][0] < this->xl + this->cell_length;
                ux = position[i][0] > this->xu - this->cell_length;
                if (lx) {
                    this->forward_communication_data[0].push_back(idpos); 
                } else if (ux) {
                    this->forward_communication_data[1].push_back(idpos);
                }
            }
        }

        void pack_ghost_particles_Y(){
            // Pack this subdomain's ghost particles into messages to its neighbouring subdomains in the Y direction
            // Here we iterate over all particles in the subdomain, including the ghost particles accumulated during the X - step
            // Accumulate the ghost particles in the communication arrays
            int3double idpos;
            bool ly, uy;
            for (int i = 0; i < this->no_particles; i++) {
                // Check for ghost particles in each direction and pack them accordingly. 
                // Particles that have moved into another subdomain are communicated here   

                idpos = int3double(i, this->positions[i]);
                // Check if the particle is across y-boundaries - these are not bounded from one direction, then *including previous ghost particles*.
                ly = position[i][0] < this->yl + this->cell_length;
                uy = position[i][0] > this->yu - this->cell_length;
                if (ly) {
                    this->forward_communication_data[2].push_back(idpos); 
                } else if (uy) {
                    this->forward_communication_data[3].push_back(idpos);
                }
            }
        }

        void pack_ghost_particles_Z(){
            // Pack this subdomain's ghost particles into messages to its neighbouring subdomains in the Z direction
            // Here we iterate over all particles in the subdomain, including the ghost particles accumulated during the XY - steps
            // Accumulate the ghost particles in the communication arrays
            int3double idpos;
            bool lz, uz;
            for (int i = 0; i < this->no_particles; i++) {
                // Check for ghost particles in each direction and pack them accordingly. 
                // Particles that have moved into another subdomain are communicated here   

                idpos = int3double(i, this->positions[i]);
                // Check if the particle is across z-boundaries - these are not bounded from one direction, then *including previous ghost particles*.
                lz = position[i][0] < this->zl + this->cell_length;
                uz = position[i][0] > this->zu - this->cell_length;
                if (lz) {
                    this->forward_communication_data[4].push_back(idpos); 
                } else if (uz) {
                    this->forward_communication_data[5].push_back(idpos);
                }        
                
            }
        }
        
        void pack_ghost_particles() {
            // Pack this subdomain's ghost particles into messages to its neighbouring subdomains
            // Check if particles are in the correct part of the positions array
            // Clear the communication arrays
            for (int i = 0; i < 6; i++) {
                this->forward_communication_data[i].clear();
            }

            // Accumulate the ghost particles in the communication arrays
            int3double idpos;
            bool lx, ux, ly, uy, lz, uz;
            for (int i = 0; i < this->no_local_particles; i++) {
                // Check for ghost particles in each direction and pack them accordingly. 
                // Particles that have moved into another subdomain are communicated here   

                idpos = int3double(i, {this->position[i][0], this->position[i][1], this->position[i][2]});

                // Check if the particle is near the lower x-boundary
                lx = position[i][0] < this->xl + this->cell_length;
                if (lx) {
                    // Check if the particle has exited the subdomain
                    lx = position[i][0] < this->xl;
                    if (lx) {
                        // Queue the removal of the particle
                        this->removal_queue.push_back(i);
                    }
                    this->forward_communication_data[0].push_back(idpos); // Index 0 needs to be changed accordingly
                }

                // Check if the particle is near the upper x-boundary
                ux = position[i][0] > this->xu - this->cell_length;
                if (ux) {
                    // Check if the particle has exited the subdomain
                    ux = position[i][0] > this->xu;
                    if (ux) {
                        // Queue the removal of the particle
                        this->removal_queue.push_back(i);
                    }
                    this->forward_communication_data[1].push_back(idpos); // Index 1 needs to be changed accordingly
                }

                // Check if the particle is near the lower y-boundary
                ly = position[i][1] < this->yl + this->cell_length;
                if (ly) {
                    // Check if the particle has exited the subdomain
                    ly = position[i][1] < this->yl;
                    if (ly) {
                        // Queue the removal of the particle
                        this->removal_queue.push_back(i);
                    }
                    this->forward_communication_data[2].push_back(idpos); // Index 2 needs to be changed accordingly
                }

                // Check if the particle is near the upper y-boundary
                uy = position[i][1] > this->yu - this->cell_length;
                if (uy) {
                    // Check if the particle has exited the subdomain
                    uy = position[i][1] > this->yu;
                    if (uy) {
                        // Queue the removal of the particle
                        this->removal_queue.push_back(i);
                    }
                    this->forward_communication_data[3].push_back(idpos); // Index 3 needs to be changed accordingly
                }

                // Check if the particle is near the lower z-boundary
                lz = position[i][2] < this->zl + this->cell_length;
                if (lz) {
                    // Check if the particle has exited the subdomain
                    lz = position[i][2] < this->zl;
                    if (lz) {
                        // Queue the removal of the particle
                        this->removal_queue.push_back(i);
                    }
                    this->forward_communication_data[4].push_back(idpos); // Index 4 needs to be changed accordingly
                }

                // Check if the particle is near the upper z-boundary
                uz = position[i][2] > this->zu - this->cell_length;
                if (uz) {
                    // Check if the particle has exited the subdomain
                    uz = position[i][2] > this->zu;
                    if (uz) {
                        // Queue the removal of the particle
                        this->removal_queue.push_back(i);
                    }
                    this->forward_communication_data[5].push_back(idpos); // Index 5 needs to be changed accordingly
                }
            }
        }

        void communicate_ghost_particles(){
            // Pack the ghost particles into message buffers for the neighbouring subdomains in order X, Y, Z
            pack_ghost_particles_X();
            communicate_ghost_particles_X();
            unpack_ghost_particles_X();

            pack_ghost_particles_Y();
            communicate_ghost_particles_Y();
            unpack_ghost_particles_Y();

            pack_ghost_particles_Z();
            communicate_ghost_particles_Z();
            unpack_ghost_particles_Z();
            
        }
        

        }
        // ------------------------------ Auxillary functions -------------------------------
        {
        bool inside_domain(double[3] pos) {
            // Test if inside domain S
            return (pos[0] > this->bounds[0] && 
                    pos[0] < this->bounds[1] && 
                    pos[1] > this->bounds[2] && 
                    pos[1] < this->bounds[3] && 
                    pos[2] > this->bounds[4] && 
                    pos[2] < this->bounds[5]);
        }
        bool inside_ghost_domain(double[3] pos) {
            // Test if inside the extended domain S'
            return (pos[0] > this->gh_bounds[0] && 
                    pos[0] < this->gh_bounds[1] && 
                    pos[1] > this->gh_bounds[2] && 
                    pos[1] < this->gh_bounds[3] && 
                    pos[2] > this->gh_bounds[4] && 
                    pos[2] < this->gh_bounds[5]);
        }
        bool is_ghost_particle(double[3] pos) {
            // Test if a particle is inside the ghost subdomain, S' \ S
            return (pos[0] > this->gh_bounds[0] && pos[0] < this->bounds[0]) || (pos[0] > this->bounds[1] && pos[0] < this->gh_bounds[1]) || 
                   (pos[1] > this->gh_bounds[2] && pos[1] < this->bounds[2]) || (pos[1] > this->bounds[3] && pos[1] < this->gh_bounds[3]) || 
                   (pos[2] > this->gh_bounds[4] && pos[2] < this->bounds[4]) || (pos[2] > this->bounds[5] && pos[2] < this->gh_bounds[5]);
        }

        void particle_tracking() {
            // Track the number of local particles and ghost particles in this subdomain
            // Track if particles need to be forward communicated
            int no_ghost_particles = 0;
            int n = this->no_subdivisions;

            // Iterate over the ghost cells and add up all the ghost particles
            // Ghost cells are at the edges of the simulation box

            this->no_particles = this->positions.size();
            this->no_ghost_particles = 0;
            this->no_local_particles = 0;
            bool inghost; 
            bool indomain; 

            for (int i = 0; i < this->no_particles; i++) {
                inghost = inside_ghost_domain(this->positions[i]);
                indomain = inside_domain(this->positions[i]);
                if (inghost) {
                    // Potentially outside the domain
                    if (max_dist > this->side_length/2.0) {
                        // Outside the domain -> communicate
                    else {
                        // Ghost particle
                        this->no_ghost_particles += 1;
                        }
                    }
                }
            }
        }

        }
}

