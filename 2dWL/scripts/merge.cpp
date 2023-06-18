#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <cmath>
#include <assert.h>
#include <complex>
#include <math.h>
#include <iomanip>




class dataset{
public:
    std::vector<double> U;
    std::vector<double> S;
    std::vector<double> N;
    std::vector<double> U_t;
    std::vector<double> S_t;
    
    double N_target;

    double avg_start;
    double avg_end;
    int n;
    void avg_ends(int n_t)
    {
        n = n_t;
        avg_start = 0;
        avg_end = 0;
        for(int i=S_t.size()-n; i<S_t.size(); i++)
        {
            avg_end+=S_t[i];
        }
        for(int i=0; i<n; i++)
        {
            avg_start+=S_t[i];
        }
        avg_start/=n;
        avg_end/=n;
    }
    void load_coord(const char* filename)
    {
        std::ifstream lp(filename);
        assert(lp.good());

        double N_temp, U_temp, S_temp;
        while(lp.good())
        {
            std::string line;
            if(std::getline(lp, line))
            {
                std::istringstream iss(line);
                iss >> N_temp >>  U_temp >> S_temp;
                if(N_temp == N_target)
                {
                    U_t.push_back(U_temp);
                    S_t.push_back(S_temp);
                }
                U.push_back(U_temp);
                S.push_back(S_temp);
                N.push_back(N_temp);
            }
        }
    }
    void shift(double shift_t)
    {
        avg_start +=shift_t;
        avg_end +=shift_t;
        
        for(int i=0; i<S.size(); i++)
        {
            S[i]+=shift_t;
        }
    }
    void dump(std::ostream& os)
    {
        for(int i=0; i<S.size()-n; i++)
            os << N[i] <<"\t" << U[i] << "\t" << std::setprecision(32) << S[i] << std::endl;
    }
    
};

int main()
{

    std::vector<std::string> filenames;
    
    filenames.push_back("../m690m600/ln_dos_1.0000000075.txt");
    filenames.push_back("../m600m500/ln_dos_1.0000000075.txt");

    
    std::vector<dataset> data;
    
    
    for(int ff=0; ff<filenames.size(); ff++)
    {
        dataset temp_data;
        temp_data.N_target = 128;

        temp_data.load_coord(filenames[ff].c_str());
        temp_data.avg_ends(5);
        data.push_back(temp_data);
    }
    
    double end = data[0].avg_end;
    double start = 0;
    for(int ff=0; ff<filenames.size(); ff++)
    {
        start = data[ff].avg_start;
        data[ff].shift(-start);
        std::cout << data[ff].avg_start << std::endl;
        data[ff].shift(end);
        end = data[ff].avg_end;
        
    }
    std::ofstream output_S("S_all.txt");
    for(int ff=0; ff<filenames.size(); ff++)
    {
        data[ff].dump(output_S);
    }
	
    return 0;
}
