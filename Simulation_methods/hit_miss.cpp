# include <cmath>
# include <iostream>
# include <random>
# include <chrono>
# include <fstream>
# define PI 3.14159565

using namespace std;

double G(double x) 
{
    double r;
    r = sqrt(1.0 - x * x);
    return r;
}

int main(void)
{
    double err, x, y;
    srand(time(NULL));

    ofstream file("Errores.txt");
    file << "M\terr" << endl;

    size_t n_in;
    auto start = chrono::high_resolution_clock::now();
    for (size_t M = 10; M < 1e9; M = 10 * M)
    {
        err = 0.0;
        for (size_t n = 0; n < 20; n++)
        {
            n_in = 0;
            for (size_t i = 0; i < M; i++)
            {
                x = (double)random() / RAND_MAX;
                y = (double)random() / RAND_MAX;

                if (y < G(x)) n_in++;
            }

            err += abs((double)n_in / M - PI / 4);
        }
        
        file << M << "\t" << err/20.0 << endl;
        cout << "M = " << M << " done" << endl;
    }

    file.close();
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);

    cout << "Time taken: " << duration.count() << " seconds" << endl;
    return 0;
}