# include <cmath>
# include <iostream>
# include <random>
# include <chrono>
# include <fstream>
# include <vector>
# include <functional>
# include <array>
# define PI 3.14159565

using namespace std;

double CDF_inv(double u)
{
    double result = - log(u);

    return result;
}

double f(double x)
{
    double result = sqrt(x) * cos(x);

    return result;
}
array<double, 2> importance(function<const double(double)> f, function<const double(double)> CDF_inv, size_t M)
{
    double mean, mean2;
    
    mean = 0.0;
    mean2 = 0.0;
    for (size_t i = 0; i < M; i++)
    {
        // Generamos valores obtenidos de una distribución exponencial
        double u = (double)random() / (exp2(31)-1);
        double x = CDF_inv(u);

        double fx = f(x);
        mean += fx;
        mean2 += fx * fx;
    }

    mean /= M;
    mean2  /= M;

    double sigma = sqrt((mean2 - mean * mean) / M);

    return {mean, sigma};
}

int main(void)
{
    int M;

    srand(time(NULL));

    cout << "Type the number of samples (M > 0): ";
    cin >> M;

    while (M < 1)
    {
        cout << "Incorrect format!" << endl;
        cout << "Type the number of samples (M > 0): ";
        cin >> M;
    }

    auto start = chrono::high_resolution_clock::now();
    auto [mean, sigma] = importance(f, CDF_inv, M);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);

    cout << "Importance sampling method: " << mean << " \u00B1 " << sigma << ", " << duration.count() << " milliseconds" << endl;
}