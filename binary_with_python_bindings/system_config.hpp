#ifndef _system_config_hpp
#define _system_config_hpp

class configuration{
public:
    int N_particles;
    std::vector < std::vector < double> > v_xyz;
    std::vector<int> v_type;
    
    void init(int N_particles_t)
    {
        N_particles = N_particles_t;
        v_xyz.clear();
        v_type.clear();
        for(int i=0; i<N_particles; i++)
        {
            std::vector < double> x_temp;
            x_temp.push_back(0);
            x_temp.push_back(0);
            x_temp.push_back(0);
            v_xyz.push_back(x_temp);
            v_type.push_back(0);
        }
        
    }
    void print_config()
    {
        for(int i=0; i<N_particles; i++)
        {
            std::cout << v_type[i] << "\t" << v_xyz[i][0] << "\t" <<  v_xyz[i][1] << "\t" << v_xyz[i][2] << std::endl;
            
        }
    }
    void perturb()
    {
        for(int i=0; i<N_particles; i++)
        {
            v_xyz[i][0] +=1.0;
        }
        
    }
};


extern void init_configuration(int N_particles_t);
//extern void set_configuration(std::vector < std::vector < double> > v_xyz_t, std::vector<int> v_type_t);
void set_configuration(coordlist_t v_xyz_t, std::vector<int> v_type_t);
extern void get_configuration(std::vector < std::vector < double> >& v_xyz_t, std::vector<int>& v_type_t);

#endif