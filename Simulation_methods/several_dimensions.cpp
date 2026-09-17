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

double cos2(const vector<double>& x)
{
    double sum2 = 0;
    double result;
    size_t l;

    l = size(x);
    for (size_t i = 0; i < l; i++)
    {
        sum2 += x[i] * x[(i+1)%l];
    }

    result = cos(sum2) * cos(sum2);

    return result;
}

double f(const vector<double>& x)
{
    double sum1 = 0;
    double sum2 = 0;
    double result;
    size_t l;

    l = size(x);
    for (size_t i = 0; i < l; i++)
    {
        sum1 -= x[i] * x[i];
        sum2 += x[i] * x[(i+1)%l];
    }

    result = exp(sum1) * cos(sum2) * cos(sum2);

    return result;
}

array<double, 2> hit_miss(function<double(const vector<double>&)> f , double xmin, double xmax, double ymax, size_t M, size_t dim)
{
    size_t n_in = 0;
    vector<double> x(dim);
    double y, I, s, mult, p;

    for (size_t i = 0; i < M; i++)
    {
        for (size_t j = 0; j < dim; j++)
        {
            x[j] = (xmax - xmin) * (double)random()/RAND_MAX + xmin;
        }
        y = ymax * (double)random()/RAND_MAX;
        
        if (y < f(x))
        {
            n_in++;
        }
        
    }

    mult = 1.0;
    for (size_t i = 0; i < dim; i++)
    {
        mult *= xmax - xmin;
    }
    
    p = (double)n_in / M;
    I = p * mult * ymax;
    s =  mult * ymax * sqrt(p * (1-p) / M);

    return {I, s};
}

array<double, 2> uniform_sampling(function<double(const vector<double>&)> f, double xmin, double xmax, size_t M, size_t dim)
{
    double f_mean, f_mean2, I, s, mult;
    vector<double> x(dim);

    f_mean = 0.0;
    f_mean2 = 0.0;
    for (size_t i = 0; i < M; i++)
    {
        for (size_t j = 0; j < dim; j++)
        {
            x[j] = (xmax - xmin) * (double)random()/RAND_MAX + xmin;
        }

        double fx = f(x);
        f_mean += fx;
        f_mean2 += fx * fx;
    }

    f_mean = f_mean / M;
    f_mean2 = f_mean2 / M;
    mult = 1.0;
    for (size_t i = 0; i < dim; i++)
    {
        mult *= xmax - xmin;
    }
    

    double var = f_mean2 - f_mean * f_mean;
    I = f_mean * mult;
    s = mult * sqrt(var / M);

    return {I, s};
}

int main(void)
{
    double xmin, xmax, ymax;
    int M, dim;

    srand(time(NULL));
    xmin = -1.0;
    xmax = 1.0;
    ymax = 1.0;
    M = 10000;

    cout << "Type the dimension (dim > 0): ";
    cin >> dim;

    while (dim < 1)
    {
        cout << "Incorrect format!" << endl;
        cout << "Type the dimension (dim > 0): ";
        cin >> dim;
    }

    cout << "Type the number of samples (M > 0): ";
    cin >> M;

    while (dim < 1)
    {
        cout << "Incorrect format!" << endl;
        cout << "Type the number of samples (M1 > 0): ";
        cin >> dim;
    }
    

    auto start1 = chrono::high_resolution_clock::now();
    auto [m_hit_miss, s_hit_miss] = hit_miss(f, xmin, xmax, ymax, M, dim);
    auto stop1 = chrono::high_resolution_clock::now();
    auto duration1 = chrono::duration_cast<chrono::milliseconds>(stop1 - start1);

    auto start2 = chrono::high_resolution_clock::now();
    auto [m_uniform, s_uniform] = uniform_sampling(f, xmin, xmax, M, dim);
    auto stop2 = chrono::high_resolution_clock::now();
    auto duration2 = chrono::duration_cast<chrono::milliseconds>(stop2 - start2);

    cout << "Hit-miss method: " << m_hit_miss << " \u00B1 " << s_hit_miss << ", " << duration1.count() << " milliseconds" << endl;
    cout << "Uniform sampling method: " << m_uniform << " \u00B1 " << s_uniform << ", " << duration2.count() << " milliseconds" << endl;
}