# include <iostream>
# include <vector>
# include <cmath>
# include <random>
# include <fstream>
using namespace std;

double magnetization(vector<vector<int>> &grid, int L)
{
    double M = 0.0;
    for (size_t i = 0; i < L; i++)
    {
        for (size_t j = 0; j < L; j++)
        {
            M += grid[i][j];
        }
    }
    
    return M / (L * L);
}

double energy(vector<vector<int>> &grid, int L)
{
    double E = 0.0;
    for (size_t i = 0; i < L; i++)
    {
        for (size_t j = 0; j < L; j++)
        {
            E -= grid[i][j] * (grid[i][(j+1)%L] + grid[(i+1)%L][j]);
        }
    }
    
    return E / (L * L);
}

void markov_step(vector<vector<int>> &grid, double T, uniform_real_distribution<double> &uni, uniform_int_distribution<int> &pos, mt19937 &rng, int L)
{
    int x = pos(rng);
    int y = pos(rng);
    int neighboors = grid[x][(y+1)%L] + grid[(x+1)%L][y] + grid[(x-1+L)%L][y] + grid[x][(y-1+L)%L];
    double p = exp(-2.0/T * grid[x][y] * neighboors);

    if (uni(rng) < p) 
    {
        grid[x][y] = -grid[x][y];
    }

    return;
}

int main(void)
{
    // Creamos el generador de números aleatorios
    // Seleccionamos la temperatura, el tamaño del sistema y el número de pasos Montecarlo
    double T = 4;
    int L = 32;
    int MCSteps = 1000;
    double E = 0.0;
    double E2 = 0.0;
    double M = 0.0;
    double M2 = 0.0;

    // Creamos los archivos de salida
    ofstream MFile("out/Magnetization.txt");
    ofstream EFile("out/Energy.txt");

    // Creamos el generador y las distribuciones
    mt19937 rng(random_device{}());
    uniform_real_distribution<double> uni(0.0, 1.0);
    uniform_int_distribution<int> pos(0, L-1);
    bernoulli_distribution bern(0.5);

    // Creamos un grid con todo qs y lo rellenamos de 1s y -1s aleatoriamente
    vector<vector<int>> grid(L, vector<int>(L, 1));
    for (size_t i = 0; i < L; i++)
    {
        for (size_t j = 0; j < L; j++)
        {
            grid[i][j] = bern(rng) ? 1:-1;
        }
        
    }

    // Ejecutamos MCSteps * L^2 iteraciones
    // Tras cada paso Montecarlo escribimos la magnetización y la energía del sistema
    int steps = 0;
    for (size_t i = 0; i < MCSteps; i++)
    {
        for (size_t j = 0; j < L*L; j++)
        {
            markov_step(grid, T, uni, pos, rng, L);
        }

        // Guardamos los datos de magnetización y energía por espín
        double curr_E = energy(grid, L);
        double curr_M = magnetization(grid, L);

        EFile << curr_E << endl;
        MFile << curr_M << endl;

        if (i > 100)
        {
            E += curr_E;
            E2 += curr_E * curr_E;
            M += curr_M;
            M2 += curr_M * curr_M;
            steps++;
        }
    }

    double EMean = E / steps;
    double E2Mean = E2 / steps;    
    double MMean = M / steps;
    double M2Mean = M2 / steps;

    double C_H = L * L * (E2Mean - EMean * EMean) / (T * T);
    double chi = L * L * (M2Mean - MMean * MMean) / T;

    EFile << C_H;
    MFile << chi;

    EFile.close();
    MFile.close();
    return 0;
}