# include <iostream>
# include <random>
# include <vector>
# include <array>
# include <fstream>

using namespace std;
struct Parameters
{
    double a;
    double b;
    double D;
    double x_0;
    double h;
    double total_time;
    double save_h;
    int N;
};

void trajectories(Parameters params)
{
    int save_each = (int)(params.save_h / params.h);
    int total_steps = (int)(params.total_time / params.h);
    int saved_steps = (int)(params.total_time / params.save_h);

    vector<double> x(params.N, params.x_0);
    ofstream file;

    // Statistical parameters
    vector<double> x_mean(saved_steps, params.x_0);
    vector<double> x2_mean(saved_steps, params.x_0 * params.x_0);
    vector<double> x4_mean(saved_steps, params.x_0 * params.x_0 * params.x_0 * params.x_0);
    vector<double> x2_variance(saved_steps, 0.0);

    file.open("Trajectories.txt");
    if (!file.is_open())
    {
        cerr << "No se pudo abrir el archivo Trajectories.txt" <<  endl;
        return;
    }

    // RNG engines
    mt19937 rng;
    normal_distribution<double> xi;

    // 1.- TRAJECTORIES
    // Escribimos el primer punto de la trayectoria
    for (size_t i = 0; i < params.N; i++)
    {
        file << x[i] << "\t";
    }
    file << x_mean[0] << "\t";
    file << x2_mean[0] << "\t";
    file << x4_mean[0] << "\t";
    file << x2_variance[0] << "\t";
    file <<  endl;
    
    // Recorremos todos los pasos que se guardan
    for (size_t t = 1; t < saved_steps; t++)
    {
        // Cada save_each pasos guardamos los datos
        for (size_t inner_t = 0; inner_t < save_each; inner_t++)
        {
            for (size_t i = 0; i < params.N; i++)
            {
                x[i] = x[i] + params.h * x[i] * (params.a - params.b * x[i] * x[i]) + sqrt(params.D * params.h) * xi(rng);
                if (inner_t == 0)
                {
                    file << x[i] << "\t";
                    x_mean[t] += x[i];
                    x2_mean[t] += x[i] * x[i];
                    x4_mean[t] += x[i] * x[i] * x[i] * x[i];
                }
            }
        }

        // Dividimos entre el número de trayectorias y guardamos en el archivo
        x_mean[t] = x_mean[t] / params.N;
        x2_mean[t] = x2_mean[t] / params.N;
        x4_mean[t] = x4_mean[t] / params.N;
        x2_variance[t] = x4_mean[t] - x2_mean[t] * x2_mean[t];
        file << x_mean[t] << "\t";
        file << x2_mean[t] << "\t";
        file << x4_mean[t] << "\t";
        file << x2_variance[t] << "\t";
        file << endl;
    }

    file.close();
}

double mean_first_passage_time(Parameters params, double R)
{
    double x, pre_x;
    double mean_fpt;
    int step = 0;

    // RNG engines
    mt19937 rng;
    normal_distribution<double> xi;
    
    // Recorremos todos los pasos que se guardan
    mean_fpt = 0.0;
    for (size_t i = 0; i < params.N; i++)
    {
        step = 0;
        pre_x = params.x_0;
        x = params.x_0;
        while (abs(x) < R)
        {

            pre_x = x;
            x = x + params.h * x * (params.a - params.b * x * x) + sqrt(params.D * params.h) * xi(rng);
            step++;
        }

        mean_fpt += params.h * (step - 1) + params.h * (R-pre_x)/(x-pre_x);
    }
    
    mean_fpt = mean_fpt / params.N;

    return mean_fpt;
}

int main()
{
    /*
    SDE: dx/dt = a * x - b * x^3 + sqrt(D) * xi(t)
    */
   
    // System parameters
    double a = 4.0;
    double b = 1.0;
    double D = 0.01; // Difussion coefficient
    double x_0 = 0.0; // Initial value
    int N = 20000;

    // Simulation parameters and variables
    double total_time = 4.0;
    double save_h = 0.1;
    double h = 0.001; // Time step
    Parameters params = {a, b, D, x_0, h, total_time, save_h, N};
    vector<double> D_list = {1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9, 1.0e-10};
    vector<double> mean_fpt(D_list.size());
    double R = 0.5;
    ofstream mfpsFile;
    mfpsFile.open("Mean_passage_time.txt");

    for (int i = 0; i < 10; i++)
    {
        params.D = D_list[i];
        mean_fpt[i] = mean_first_passage_time(params, R);
        cout << "D: " << D_list[i] << " | mfpt: " << mean_fpt[i] << endl;
        mfpsFile << D_list[i] << "\t" << mean_fpt[i] << endl;
    }

    mfpsFile.close();
    return 0;
}