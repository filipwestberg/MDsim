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

struct IDPOSITION
{
    int id;
    double position[3];
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

