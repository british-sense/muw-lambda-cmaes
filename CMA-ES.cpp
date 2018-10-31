#include <random>
#include <climits>
#include "algebra.hpp"

using namespace std;
using array = vector<double>;
using matrix = vector<vector<double> >;

random_device seed;
mt19937 mt(seed()); 
normal_distribution<> dist(0., 1.);

static const int n = 100;
static const double upper_limit = 5.12;
static const double lower_limit = -5.12;

class individual{
    
    public:
    array x;
    array z;
    double fitness;

    public:
    individual() : x(zeros(n)), z(zeros(n)), fitness(0.){
    }
    void evaluate(){
        // rastrigin function
        fitness = 10. * n;
        for(int i = 0; i < n; i++) fitness += x[i] * x[i] - 10. * cos(2. * M_PI * x[i]);
    }
    bool operator<(const individual& rhs)const{
        return fitness < rhs.fitness;
    }
    bool operator>(const individual& rhs)const{
        return fitness > rhs.fitness;
    }
};

int main(){

    array x_mean_weight(zeros(n));
    uniform_real_distribution<> range(lower_limit, upper_limit);
    for(int i = 0; i < n; i++) x_mean_weight[i] = range(mt);

    array z_mean_weight(n);
    double sigma = 100.;
    double min_sigma = 1e-15;

    const double lambda = 4. + floor(3. * log(n));
    const double mu = floor(lambda / 2.);

    array array_weight(ones(n));// = log((lambda + 1) / 2) - log(seq(1, mu));

    double c_c = 4. / (n + 4.);
    double c_cov = 2. / ((n + sqrt(2.)) * (n + sqrt(2.)));
    double c_s = 4. / (n + 4.);
    double damp = (1. / c_s) + 1.;

    matrix B(identity(n));
    matrix D(identity(n));
    matrix BD = B * D;
    matrix C = BD * trans(BD);
    array p_c(zeros(n));
    array p_s(zeros(n));
    double c_w = sum(array_weight) / norm(array_weight);
    const double chi_n = sqrt(n) * (1. - (1. / (4. * n)) + (1. / (21. * n * n)));

    int count_evaluate = 0;
    int generation = 1000;
    double min = INT_MAX;

    for(int t = 1; t <= generation; t++){

        // 個体の生成
        vector<individual> population(lambda);
        for(int k = 0; k < lambda; k++){
            for(int i = 0; i < n; i++) population[k].z[i] = dist(mt);
            population[k].x = x_mean_weight + sigma * (BD * population[k].z); 
            population[k].evaluate();
            count_evaluate++;
        }

        // 評価値でソート
        sort(population.begin(), population.end());

        // 平均ベクトルの更新
        x_mean_weight = zeros(x_mean_weight.size());
        z_mean_weight = zeros(z_mean_weight.size());
        for(int k = 0; k < mu; k++){
            x_mean_weight += population[k].x * array_weight[k];
            z_mean_weight += population[k].z * array_weight[k];
        }
        x_mean_weight /= sum(array_weight);
        z_mean_weight /= sum(array_weight);


        if(min > population.front().fitness) min = population.front().fitness;
        cout << t << " : " << min << endl;

        // パスの更新
        p_c = (1. - c_c) * p_c + (sqrt(c_c * (2. - c_c)) * c_w) * (BD * z_mean_weight);
        C = (1. - c_cov) * C + (c_cov * p_c * p_c);
        p_s = (1. - c_s) * p_s + (sqrt(c_s * (2. - c_s)) * c_w) * (B * z_mean_weight);
        sigma = sigma * exp(((norm(p_s) - chi_n) / chi_n) / damp);
        
        if(((int)(count_evaluate / lambda) % (int)(n / 10)) < 1.){
            C = upper_triangle(C) + trans(upper_triangle(C, 1));
            // diagonalization(C, B, D);
            jacobi(C, B, D);
            if(D.back().back() > 1e14 * D.front().front()){
                long long tmp = (D.back().back() / 1e14) - D.front().front();
                C += tmp * identity(n);
                D += tmp * identity(n);
            }
            std::vector<double> dv = sqrt(diagonalization_component(D));
            D = diagonalization_matrix(dv);
            BD = B * D;
        }

        // ステップサイズの調整
        if(sigma * D.front().front() < min_sigma || population.front().fitness == population[mu].fitness || x_mean_weight == (x_mean_weight + 0.2 * sigma * col(BD, 1 + floor((int)(count_evaluate / lambda) % n)))) sigma *= 1.4;
    }
    return 0;
}