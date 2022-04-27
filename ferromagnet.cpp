#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

using namespace boost::accumulators;
using namespace std;

class Heisenberg_QMC
{
public:
    int M, n;
    bool **spins;
    bool **crossbondright;
    bool **crossbondleft;
    bool **ybonds;
    int **clusters;
    double beta, J = 1.0;
    unsigned seed;
    double delta, P;
    double *Szi_Szj;
    int *properLabels;

    std::default_random_engine generator;
    std::uniform_int_distribution<int> random_i = std::uniform_int_distribution<int>(0, M - 1);
    std::uniform_int_distribution<int> random_j = std::uniform_int_distribution<int>(0, n - 1);
    std::uniform_real_distribution<double> random_P = std::uniform_real_distribution<double>(0, 1);

    accumulator_set<double, features<tag::sum>, void> Z;
    accumulator_set<double, features<tag::sum>, void> *Sij;

    void spins_to_terminal()
    {
        //prints spins to terminal

        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (j != n - 1)
                    std::cout << spins[i][j] << " ";
                else
                    std::cout << spins[i][j] << "\n";
            }
        }
        std::cout << "\n"
                  << std::flush;
    }

    Heisenberg_QMC(int M, int n, unsigned seed) : M(M), n(n), seed(seed)
    {
        spins = new bool *[M];
        ybonds = new bool *[M];
        crossbondleft = new bool *[M];
        crossbondright = new bool *[M];
        properLabels = new int[M * n];
        clusters = new int *[M];

        for (int i = 0; i < n * M; i++)
            properLabels[i] = 0;

        for (int i = 0; i < M; i++)
        {
            spins[i] = new bool[n];
            ybonds[i] = new bool[n];
            crossbondleft[i] = new bool[n];
            crossbondright[i] = new bool[n];
            clusters[i] = new int[n];
        }
        generator.seed(seed);
    }

    std::string tile_type(int i, int j)
    {
        //returns tile type
        //if (i==M) i=0;
        //if (j==n) j=0;
        //if (i==-1) i=M-1;
        //if (j==-1) j = n-1;

        int right = (j + 1) % n;
        int up = (i + 1) % M;

        if (spins[i][j] == 1 && spins[i][right] == 0 && spins[up][right] == 0 && spins[up][j] == 1)
        {
            return "c1";
        }
        else if (spins[i][j] == 0 && spins[i][right] == 1 && spins[up][right] == 1 && spins[up][j] == 0)
        {
            return "c2";
        }
        else if (spins[i][j] == 1 && spins[i][right] == 1 && spins[up][right] == 1 && spins[up][j] == 1)
        {
            return "a1";
        }
        else if (spins[i][j] == 0 && spins[i][right] == 0 && spins[up][right] == 0 && spins[up][j] == 0)
        {
            return "a2";
        }
        else if (spins[i][j] == 1 && spins[i][right] == 0 && spins[up][right] == 1 && spins[up][j] == 0)
        {
            return "b1";
        }
        else if (spins[i][j] == 0 && spins[i][right] == 1 && spins[up][right] == 0 && spins[up][j] == 1)
        {
            return "b2";
        }
        else
        {
            std::cout << i << " " << j << "\n";
            return "not valid";
        }
    }

    double weight(std::string type)
    {
        //returns weight of a tile type
        if (type == "a1" || type == "a2")
        {
            return exp(+delta / 4) * 0.9;
        }
        else if (type == "b1" || type == "b2")
        {
            return exp(-delta / 4) * cosh(delta / 2) * 0.9;
        }
        else if (type == "c1" || type == "c2")
        {
            return exp(-delta / 4) * sinh(delta / 2) * 0.9;
        }
        else
        {
            //spins_to_terminal();
            //std::cin.get();
            return 0;
            std::cout << "problem weight";
        }
    }
    void make_ybond(int i, int j)
    {
        //connects spins in y direction
        ybonds[i][j] = true;
        ybonds[i][(j + 1) % n] = true;
    }
    void make_crossbond(int i, int j)
    {
        //connects diagonal spins
        crossbondright[i][j] = true;
        crossbondleft[i][(j + 1) % n] = true;
    }
    void make_bonds()
    {
        for (int i = 0; i < M; i++)
        {
            for (int j = (i % 2); j < n; j = j + 2)
            {
                std::string type = tile_type(i, j);
                if (type == "not valid")
                    std::cout << "no bonds\n";
                if (type == "a1" || type == "a2")
                {
                    if (random_P(generator) < weight("c1") / (weight("c1") + weight("b1")))
                    {
                        //if (random_P(generator)<(0)){
                        //if (random_P(generator)<0.5*delta/(1+0.25*delta)){
                        make_crossbond(i, j);
                    }
                    else
                        make_ybond(i, j);
                }
                else if (type == "c1" || type == "c2")
                {
                    make_ybond(i, j);
                }
                else if (type == "b1" || type == "b2")
                {
                    make_crossbond(i, j);
                }
            }
        }
    }
    int properLabel(int label)
    {
        //helper function for labeling loops
        while (properLabels[label] != label)
        {
            label = properLabels[label];
        }
        return label;
    }

    double calc_Z()
    {
        //calculates statistical weight of a state
        double Z = 1.0;
        for (int i = 0; i < M; i++)
        {
            for (int j = i % 2; j < n; j = j + 2)
            {
                Z = Z * weight(tile_type(i, j));
            }
        }
        return Z;
    }

    void label_clusters()
    {
        //does all the labeling
        int label = 0;
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < n; j++)
            {
                int ibond = -7000;
                int jbond = -7001;
                int bonds = 0;
                if (ybonds[i][j] == true)
                {
                    ibond = (i + 1) % M;
                    jbond = j;
                    bonds++;
                }
                if (crossbondright[i][j] == true)
                {
                    ibond = (i + 1) % M;
                    jbond = (j + 1) % n;
                    bonds++;
                }
                if (crossbondleft[i][j] == true)
                {
                    ibond = (i + 1) % M;
                    jbond = (j - 1);
                    if (jbond == -1)
                        jbond = n - 1;
                    bonds++;
                }

                if (clusters[i][j] == -1 && i != M - 1)
                {
                    clusters[i][j] = label;
                    properLabels[label] = label;
                    label++;
                }
                if (i != M - 1)
                    clusters[ibond][jbond] = clusters[i][j];
            }
        }

        for (int k = 0; k < n; k++)
        {
            int i = M - 1;

            for (int j = 0; j < n; j++)
            {
                int bonds = 0;
                int ibond, jbond;
                if (ybonds[i][j] == true)
                {
                    ibond = (i + 1) % M;
                    jbond = j;
                    bonds++;
                }
                if (crossbondright[i][j] == true)
                {
                    ibond = (i + 1) % M;
                    jbond = (j + 1) % n;
                    bonds++;
                }
                if (crossbondleft[i][j] == true)
                {
                    ibond = (i + 1) % M;
                    jbond = (j - 1);
                    if (jbond == -1)
                        jbond = n - 1;
                    bonds++;
                }

                if (bonds != 1)
                    std::cout << "problemmm2";
                int plabel1 = properLabel(clusters[i][j]);
                int plabel2 = properLabel(clusters[0][jbond]);
                if (plabel1 <= plabel2)
                    properLabels[clusters[i][j]] = plabel1;
                else
                    properLabels[clusters[i][j]] = plabel2;
            }
        }
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < n; j++)
            {
                clusters[i][j] = properLabel(clusters[i][j]);
            }
        }
    }

    void flip_clusters()
    {
        //randomly flips clusters
        int ncluster = 0;
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (clusters[i][j] > ncluster)
                    ncluster = clusters[i][j];
            }
        }
        for (int m = 0; m < ncluster + 1; m++)
        {
            bool flip = false;
            if (random_P(generator) > 0.5)
            {
                flip = true;
            }
            if (flip)
            {
                for (int i = 0; i < M; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        if (clusters[i][j] == m)
                            spins[i][j] = !spins[i][j];
                    }
                }
            }
        }
    }

    double calcZ()
    {
        double Z = 1;
        for (int i = 0; i < M; i++)
        {
            for (int j = (i % 2); j < n; j = j + 2)
            {
                Z = Z * weight(tile_type(i, j));
            }
        }
        return Z;
    }

    void update_Z()
    {
        Z(calcZ());
    }

    void step()
    {
        //performs one monte carlo step
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < n; j++)
            {
                clusters[i][j] = -1;
                ybonds[i][j] = false;
                crossbondright[i][j] = false;
                crossbondleft[i][j] = false;
            }
        }
        for (int i = 0; i < M * n; i++)
        {
            properLabels[i] = 0;
        }

        make_bonds();

        label_clusters();
        flip_clusters();
    }

    void thermalize_spins()
    {
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < n; j++)
            {
                spins[i][j] = false;
                clusters[i][j] = -1;
            }
        }
        for (int k = 0; k < M * n; k++)
        {
            step();
            //spins_to_terminal();
            //std::cin.get();
        }
    }

    void set_params(double beta)
    {
        delta = 2.0 * beta / M;
    }

    void runSim(double beta, int steps)
    {
        //runs the simulation
        set_params(beta);
        thermalize_spins();
        Z = {};

        for (int i = 0; i < steps; i++)
        {
            step();
            update_Z();
        }
    }

    void run_Sim_output(int steps, std::string filename)
    {
        //runs the sim and outputs data
        double beta_start = 0.2;
        double beta_end = 10;
        int beta_steps = 50;
        std::ofstream out;
        out.open(filename);
        for (int i = 0; i < beta_steps; i++)
        {
            double beta = beta_start + (beta_end - beta_start) / beta_steps * i;
            std::cout << "beta: " << beta << "\n";
            runSim(beta, steps);
            out << beta << " " << sum(Z) << "\n";
        }
        out.close();
    }
};

void TemperatureDependence_ensemble(std::string filename, int M, int n, int steps, double bstart, double bend, int num_states)
{
    //runs multiple simulations using openmp
    int bsteps = 100;
    std::ofstream out;
    out.open(filename);

    Heisenberg_QMC **p = new Heisenberg_QMC *[num_states];
    for (int i = 0; i < num_states; i++)
    {
        p[i] = new Heisenberg_QMC(M, n, time(0) + 2754 * i);
    }

    for (int i = 0; i < bsteps + 1; i++)
    {

        double beta = bstart + i * (bend - bstart) / bsteps;
        std::cout << "beta=" << beta << "\n";

#pragma omp parallel for
        for (int j = 0; j < num_states; j++)
        {
            p[j]->runSim(beta, steps);
            std::cout << "sim " << j << " finished\n";
        }
        std::cout << "output!\n";
        double Eavg = 0.0;
        for (int k = 0; k < num_states; k++)
        {
            Eavg += sum(p[k]->Z);
        }
        Eavg /= num_states;

        out << beta << " " << Eavg << "\n";
    }
    out.close();
}

int main(void)
{

    //example for n=12 spin chain
    int M = 10;
    int n = 12;
    int steps = 2000000;
    double Tstart = 0.1;
    double Tend = 10;
    int num_states = 4;
    TemperatureDependence_ensemble("output12test.txt", M, n, steps, 1.0 / Tend, 1.0 / Tstart, num_states);
    //Heisenberg_QMC hm = Heisenberg_QMC(M, n, time(0)+2754);
    //hm.run_Sim_output(200000, "testdata2.txt");
    return 0;
}
