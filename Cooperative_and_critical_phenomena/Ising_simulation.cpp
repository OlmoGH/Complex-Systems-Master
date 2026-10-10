# include <iostream>
# include <vector>
# include <cmath>
# include <random>
# include <fstream>
# include <string>
# include <sstream>
# include <iomanip>
# include <array>

using namespace std;

struct Solution {
    double T;
    double E;
    double M;
    double C_H;
    double chi;
};

inline int get_id(const int x, const int y, const int L, const int mask)
{
    return (x & mask) + (y & mask) * L;
}

vector<int> initialise(const int L, bernoulli_distribution &bern, mt19937 &rng)
{
    vector<int> grid(L * L, 1);
    for (int i = 0; i < L * L; i++)
    {
        grid[i] = bern(rng) ? 1:-1;
    }

    return grid;
}

double magnetization(vector<int> &grid, int L, double inv_N)
{
    double M = 0.0;
    for (int i = 0; i < L * L; i++)
    {
        M += grid[i];

    }
    
    return M * inv_N;
}

double energy(const vector<int> &grid, const int L, double inv_N, const int mask)
{
    double E = 0.0;
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            E -= grid[get_id(i, j, L, mask)] * (grid[get_id(i, j+1, L, mask)] + grid[get_id(i+1, j, L, mask)]);
        }
    }
    
    return E * inv_N;
}

void metropolis_step(
    vector<int> &grid, 
    const double T, 
    double &curr_E, 
    double &curr_M, 
    uniform_real_distribution<double> &uni, 
    uniform_int_distribution<int> &pos, 
    mt19937 &rng, 
    const int L, 
    const double prob8, 
    const double prob4,
    const double inv_N,
    const int mask)
{
    int x = pos(rng);
    int y = pos(rng);
    int id = get_id(x, y, L, mask);
    int right = grid[get_id(x, y+1, L, mask)];
    int down = grid[get_id(x+1, y, L, mask)];
    int up = grid[get_id(x-1, y, L, mask)];
    int left = grid[get_id(x, y-1, L, mask)];
    int neighbors = right + down + up + left;
    int dE = 2 * grid[id] * neighbors;
    bool accept = true;

    if (dE > 0) {
        // Si el cambio de energía es 4 comprobamos si p < exp(-4/T) y si es 8 comprobamos si p < exp(-8/T)
        accept = dE == 4 ? uni(rng) < prob4 : uni(rng) < prob8;
    }

    if (accept)
    {
        grid[id] = -grid[id];
        curr_E += (double)dE * inv_N;
        curr_M += 2 * (double)grid[id] * inv_N;
    }

    return;
}

int wolff_step(
    vector<int> &grid, 
    double &curr_E, 
    double &curr_M, 
    uniform_real_distribution<double> &uni, 
    uniform_int_distribution<int> &pos, 
    mt19937 &rng, 
    const int L, 
    const double prob, 
    vector<int> &stack,
    const double inv_N,
    const int mask)
{
    int x = pos(rng);
    int y = pos(rng);
    int sign = grid[get_id(x, y, L, mask)];
    int cluster_size = 1;

    // Guardamos el valor del espín inicial y colocamos el espín inicial en la cola
    grid[get_id(x, y, L, mask)]  = -sign;
    stack.push_back(get_id(x, y, L, mask));

    // Recorremos los espines elegibles
    while(!stack.empty())
    {
        auto actual_pos = stack.back();
        int x_actual = actual_pos % L;
        int y_actual = actual_pos / L;
        stack.pop_back();

        // Seleccionamos los cuatro vecinos y comprobamos si son elegibles
        int neighbors[4] = {
            get_id(x_actual + 1, y_actual, L, mask), 
            get_id(x_actual - 1, y_actual, L, mask), 
            get_id(x_actual, y_actual + 1, L, mask), 
            get_id(x_actual, y_actual - 1, L, mask)
        };

        for (const auto &actual_pos : neighbors)
        {
            if (grid[actual_pos] == sign)
            {
                if (uni(rng) < prob)
                {
                    grid[actual_pos] = -sign;
                    stack.push_back(actual_pos);
                    cluster_size++;
                }
            }
        }
    }

    curr_M += -2.0 * cluster_size * sign * inv_N;
    
    return cluster_size;
}

inline void move_spin(
    const int pos, 
    const int alligned,
    array<vector<int>, 5> &spin_types, 
    vector<int> &spin_class, 
    vector<int> &spin_pos)
{
    int old_class = spin_class[pos];
    int new_class = old_class + alligned;
    if (new_class == old_class) return;

    int old_pos = spin_pos[pos];
    auto aux = spin_types[old_class].back(); // Seleccionamos el últmo espín de la clase del nuestro
    spin_types[old_class][old_pos] = aux; // Movemos el último espín a la posición de nuestro espín 
    spin_pos[aux] = spin_pos[pos]; // Actualizamos la posición del espín aux
    spin_types[old_class].pop_back(); // Eliminamos el último espín de la antigua clase
    spin_types[new_class].push_back(pos);
    spin_pos[pos] = spin_types[new_class].size() - 1;
    spin_class[pos] = new_class;

    return;
}

void classificate_spins(
    vector<int> &grid, 
    array<vector<int>, 5> &spin_types, 
    vector<int> &spin_class, 
    vector<int> &spin_pos,
    const int L)
    {
        int mask = L-1;
        for (int i = 0; i < L; i++)
        {
            for (int j = 0; j < L; j++)
            {
                int energy = 0;

                energy += grid[get_id(i+1, j, L, mask)];
                energy += grid[get_id(i-1, j, L, mask)];
                energy += grid[get_id(i, j+1, L, mask)];
                energy += grid[get_id(i, j-1, L, mask)];
                energy *= grid[get_id(i, j, L, mask)];

                int id = (4 - energy) / 2;
                int pos = get_id(i, j, L, mask);
                spin_types[id].push_back(pos);
                spin_class[pos] = id;
                spin_pos[pos] = spin_types[id].size()-1;

            }
            
        }
        
    }


double BKL_step(
    vector<int> &grid, 
    array<vector<int>, 5> &spin_types, 
    vector<int> &spin_class,
    vector<int> &spin_pos,
    double &curr_E, 
    double &curr_M,
    uniform_real_distribution<double> &uni, 
    uniform_int_distribution<int> &pos, 
    mt19937 &rng, 
    const int L, 
    const double prob4, 
    const double prob8,
    const double inv_N,
    const int mask)
{
    double acceptance[5] = {prob8, prob4, 1.0, 1.0, 1.0};
    double total_propensity = 0.0;

    for (int i = 0; i < 5; i++)
    {
        total_propensity += spin_types[i].size() * acceptance[i];
    }

    // Generamos un número aleatorio entre 0 y total_propensity
    double propensity = uni(rng) * total_propensity;

    int idx = 0;
    double sum = 0.0;
    for (; idx < 5; idx++)
    {

        sum += spin_types[idx].size() * acceptance[idx];

        if (sum >= propensity) break;
    }

    if (idx == 5) idx = 4;

    if (spin_types[idx].empty()) {
        for (idx = 0; idx < 5; idx++) {
            if (!spin_types[idx].empty()) break;
        }
    }

    // Invertimos un espín al azar del tipo seleccionado y movemos ese y sus vecinos a los slots correspondientes
    int rand_id = static_cast<int>(uni(rng) * spin_types[idx].size());
    int xyid = spin_types[idx][rand_id];
    grid[xyid] = -grid[xyid];
    int x = xyid % L;
    int y = xyid / L;

    int last_spin = spin_types[idx].back();
    spin_types[idx][rand_id] = last_spin;
    spin_pos[last_spin] = rand_id;
    spin_types[idx].pop_back();

    spin_types[4-idx].push_back(xyid);
    spin_class[xyid] = 4 - idx; // Actualizamos la clase de energía dle espín xy
    spin_pos[xyid] = spin_types[4-idx].size()-1;


    // Cambiamos la clase a la que pertenece cada espín
    int alligned_up = grid[get_id(x+1, y, L, mask)] == grid[xyid]? -1:1;
    move_spin(get_id(x+1, y, L, mask), alligned_up, spin_types, spin_class, spin_pos);

    int alligned_down = grid[get_id(x-1, y, L, mask)] == grid[xyid]? -1:1;
    move_spin(get_id(x-1, y, L, mask), alligned_down, spin_types, spin_class, spin_pos);

    int alligned_right = grid[get_id(x, y+1, L, mask)] == grid[xyid]? -1:1;
    move_spin(get_id(x, y+1, L, mask), alligned_right, spin_types, spin_class, spin_pos);

    int alligned_left = grid[get_id(x, y-1, L, mask)] == grid[xyid]? -1:1;
    move_spin(get_id(x, y-1, L, mask), alligned_left, spin_types, spin_class, spin_pos);
    
    // Actualizamos la magnetización y la energía
    curr_M += (2.0 * grid[xyid]) * inv_N;
    curr_E += (2.0 - idx) * 4.0 * inv_N;

    // Calculamos dt
    double dt = -log(1 - uni(rng)) / total_propensity;
    return dt;
}

Solution simulate_and_save(
    vector<int> &grid, 
    const double T, 
    const int L, 
    const int MCSteps, 
    const int warmup, 
    const int save_each,
    mt19937 &rng, 
    uniform_real_distribution<double> &uni, 
    uniform_int_distribution<int> &pos)
{
    // Creamos el generador de números aleatorios
    double inv_N = 1.0 / (L * L);
    const int mask = L-1;
    double E = 0.0;
    double E2 = 0.0;
    double M = 0.0;
    double M2 = 0.0;
    double curr_E = energy(grid, L, inv_N, mask);
    double curr_M = magnetization(grid, L, inv_N);
    double prob8 = exp(-8/T);
    double prob4 = exp(-4/T);
    double prob = 1.0 - exp(-2.0/T);
    double total_time = 0;
    double dt = 1.0;
    array<vector<int>, 5> spin_types;
    vector<int> spin_class(L * L);
    vector<int> spin_pos(L * L);
    vector<int> stack;
    int cluster_accumulate = 0;
    double time_accumulate = 0;
    int current_wipe = 0;
    int last_saved_wipe = warmup - save_each;
    stack.reserve(L*L);
    classificate_spins(grid, spin_types, spin_class, spin_pos, L);

    // Ejecutamos MCSteps * L^2 pasos de Metropolis
    // Tras cada paso Montecarlo escribimos la magnetización y la energía del sistema
    // Para T ~ T_crit aplicamos Wolff, para T > T_crit aplicamos metropolis y para T < T_crit aplicamos Guillespie
    
    if (T < 2.0)
    {
        while (time_accumulate < MCSteps)
        {
            dt = BKL_step(grid, spin_types, spin_class, spin_pos, curr_E, curr_M, uni, pos, rng, L, prob4, prob8, inv_N, mask);
            time_accumulate += dt;

            if (time_accumulate > warmup)
            {            
                E += curr_E * dt;
                E2 += curr_E * curr_E * dt;
                M += abs(curr_M) * dt;
                M2 += curr_M * curr_M * dt;
                total_time += dt;
            }
        }
        
        
    }
    else if (T < 3.0) 
    {
        while (cluster_accumulate < MCSteps * L * L)
        {
            stack.clear();
            cluster_accumulate += wolff_step(grid, curr_E, curr_M, uni, pos, rng, L, prob, stack, inv_N, mask);
            current_wipe = cluster_accumulate / (L * L);

            if ((current_wipe >= warmup) && (current_wipe - last_saved_wipe >= save_each))
            {
                last_saved_wipe = current_wipe;
                curr_E = energy(grid, L, inv_N, mask);
                E += curr_E * dt;
                E2 += curr_E * curr_E * dt;
                M += abs(curr_M) * dt;
                M2 += curr_M * curr_M * dt;
                total_time += 1.0;
            }
            
        }
        
    }
    else
    {
        const int total_steps = MCSteps * L * L;
        const int save_interval = save_each * L * L;
        int next_save = warmup * L * L; // Marcador del primer guardado  

        for (int i = 0; i < total_steps; i++)
        {
            metropolis_step(grid, T, curr_E, curr_M, uni, pos, rng, L, prob8, prob4, inv_N, mask);

            if (i == next_save)
            {
                E += curr_E * dt;
                E2 += curr_E * curr_E * dt;
                M += abs(curr_M) * dt;
                M2 += curr_M * curr_M * dt;
                total_time += 1.0;

                next_save += save_interval;
            }
        }
        
    }

    double EMean = E / total_time;
    double E2Mean = E2 / total_time;    
    double MMean = M / total_time;
    double M2Mean = M2 / total_time;

    double C_H = L * L * (E2Mean - EMean * EMean) / (T * T);
    double chi = L * L * (M2Mean - MMean * MMean) / T;

    // Devolvemos la temperatura, Energía, Magnetización, Capacidad calorífica y susceptibilidad magnética
    return {T, EMean, MMean, C_H, chi};
}

vector<double> correlation_function(
    vector<int> &grid, 
    const double T, 
    const int L, 
    const int MCSteps,
    mt19937 &rng, 
    uniform_real_distribution<double> &uni, 
    uniform_int_distribution<int> &pos)
{
    double inv_N = 1.0 / (L * L);
    const int mask = L-1;    
    double curr_E = energy(grid, L, inv_N, mask);
    double curr_M = magnetization(grid, L, inv_N);
    double prob8 = exp(-8/T);
    double prob4 = exp(-4/T);
    const int total_steps = MCSteps * L * L;
    vector<double> magnetization_series(MCSteps, 0.0);
    vector<double> correlation(MCSteps, 0.0);
    double mean_magnetization = 0.0;

    for (int i = 0; i < total_steps; i++)
    {
        metropolis_step(grid, T, curr_E, curr_M, uni, pos, rng, L, prob8, prob4, inv_N, mask);
        if (i % (L * L) == 0)
        {
            magnetization_series[i/(L * L)] = curr_M;
        }
    }

    // for (int i = 0; i < MCSteps; i++)
    // {
    //     mean_magnetization += magnetization_series[i] ;
    // }

    // mean_magnetization = mean_magnetization / MCSteps;
    

    // for (int i = 0; i < MCSteps; i++)
    // {
    //     for (int j = 0; j < MCSteps - i; j++)
    //     {
    //         correlation[i] += magnetization_series[j] * magnetization_series[j + i] - mean_magnetization * mean_magnetization;
    //     }

    //     correlation[i] = correlation[i] / (MCSteps - i);
    // }
    
    return magnetization_series;
}

int main(void)
{
    // Creación del archivo de salida 
    ofstream OutFile("out/Data.txt");
    ofstream CorrFile("out/Correlation.txt");
    if (!OutFile.is_open())
    {
        cout << "No se pudo abrir el archivo Data.txt" << endl;
        return -1;
    }
    if (!CorrFile.is_open())
    {
        cout << "No se pudo abrir el archivo Correlation.txt" << endl;
        return -1;
    }
    OutFile << "T\tE\tM\tC_H\tchi" << endl;

    // Creamos las variables de la simulación y el grid inicial
    double T_0 = 4.0; // Temperatura inicial del sistema
    double dT = 0.01; // Incremento de la temperatura
    double T_min = 1.0; // Temperatura mínima
    Solution data; // Vector de datos
    int L = 32; // Tamaño del sistema
    int MCSteps = 10000; // Número de pasos Montecarlo
    int warmup = 100; // Pasos de calentamiento
    int save_each = 1000; // Cada cuantos pasos se guardan los datos
    mt19937 rng(random_device{}());
    uniform_real_distribution uni(0.0, 1.0);
    uniform_int_distribution pos(0, L-1);
    bernoulli_distribution bern(0.5);
    vector<int> grid = initialise(L, bern, rng);

    auto correlation = correlation_function(grid, 4.0, L, MCSteps, rng, uni, pos);
    for (int i = 0; i < MCSteps; i++)
    {
        CorrFile << correlation[i] << endl;
    }
    
    // for (double T = T_0; T > T_min; T -= dT)
    // {
    //     cout << "Simulando sistema para T = " + to_string(T) << endl;
    //     data = simulate_and_save(grid, T, L, MCSteps, warmup, save_each, rng, uni, pos);
    //     OutFile << data.T << "\t" << data.E << "\t" << data.M << "\t" << data.C_H << "\t" << data.chi << endl;
    // }
    
    OutFile.close();
    CorrFile.close();
    return 0;
}