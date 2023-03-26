#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <utility>

using namespace std;

const int k = 10; //状態遷移数
const int n = 10; //ステップ数
int org_vx[k+1] = {  -2, -1, -1,  0,  0, 0,  1,  1,  2,  2};
int org_vy[k+1] = {  0,  0,  1, -2, -1,  0, -1,  0, -1,  0};
int sx = 1;
int sy = 1;
int cnt = 0;
double dt = 0.1;
double epsilon = 1.0;
double n_11 = 0.5;
double n_22 = 0.5;
double x_1 = 1.0;
double x_2 = 0.5;
double sum_weight = 0;

double dt_M = dt / (double)n;
void Exhaustive(int chosen[], int index, int r, int start, int end, vector<double>& weight);
void ExhaustiveRepetition(int k, int n, vector<double>& weight);

int main(int argc, char *argv[]) {

    const int r = 14;
    vector<double> weight;
    weight.resize(r+1);
    int walk_corres[r+1];
    int walk_x[r+1] = {-2, -1, -1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2};
    int walk_y[r+1] = {0, 0, 1, -2, -1, -1, -1, 0, 0, -1, -1, 0, -1, 0};
    double rate[r+1] = {0.5, 1.0, 1.0, 0.5, -1.0, -1.0, 1.0, -1.0, 1.0, -2.0, -1.0, -2.0, -1.0, -1.0};

    int index_rate1[r+1] = {2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int index_rate2[r+1] = {0, 0, 0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    int index_ep[r+1] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1};
    int index_n1[r+1] = {2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int index_n2[r+1] = {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int index_x1[r+1] = {0, 0, 0, 0, 1, 2, 0, 2, 0, 1, 0, 1, 0, 0};
    int index_x2[r+1] = {0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0};

    // for(int i=0;i<r;i++){
    //     for(int j=0;j<k;j++){
    //         if(walk_x[i]==vx[j] && walk_y[i]==vy[j]){
    //             walk_corres[i] = j;

    //         }
    //     }
    // }

    for(int i=0;i<r;i++){ //微分演算子以外の計算
        weight[i] = rate[i];
        for(int j=0;j<index_ep[i];j++){
            weight[i] *= epsilon;
        }

        for(int j=0;j<index_n1[i];j++){
            weight[i] *= n_11;
        }

        for(int j=0;j<index_n2[i];j++){
            weight[i] *= n_22;
        }

        for(int j=0;j<index_x1[i];j++){
            weight[i] *= x_1;
        }

        for(int j=0;j<index_x2[i];j++){
            weight[i] *= x_2;
        }
    }

    chrono::system_clock::time_point start, end;

    start = chrono::system_clock::now();

    ExhaustiveRepetition(k, n, weight);
    cout << std::setprecision(16) << sum_weight << endl;

    end = chrono::system_clock::now();

    // cout << "sx=" << sx << ' ' << "sy=" << sy << endl;
    // cout << "cnt=" << cnt << endl << endl;

    double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    printf("time %lf[ms]\n", time);

    return 0;
}

void Exhaustive(int chosen[], int index, int r, int start, int end, vector<double>& weight)
{
    if (index == n)
    {
        double weight_temp=1;
        double temp;
        int sum_x=sx;
        int sum_y=sy;

        for (int i = 0; i < n; i++){

            if(chosen[i] == 0){
                temp = weight[0]*sum_x*(sum_x-1);
            }else if(chosen[i] == 1){
                temp = weight[1]*sum_x;

            }else if(chosen[i] == 2){
                temp = weight[2]*sum_x;

            }else if(chosen[i] == 3){
                temp = weight[3]*sum_y*(sum_y-1);

            }else if(chosen[i] == 4){
                temp = weight[4]*sum_y;
                temp += weight[5]*sum_y;
                temp += weight[6]*sum_y;

            }else if(chosen[i] == 5){
                temp = weight[7]*sum_y;
                temp += weight[8]*sum_y;

            }else if(chosen[i] == 6){
                temp = weight[9]*sum_y;
                temp += weight[10]*sum_y;

            }else if(chosen[i] == 7){
                temp = weight[11]*sum_y;

            }else if(chosen[i] == 8){
                temp = weight[12]*sum_y;

            }else if(chosen[i] == 9){
                temp = weight[13]*sum_y;

            }

            if(chosen[i] == 5){
                weight_temp = weight_temp*(1+temp*dt_M);
            }else{
                weight_temp = weight_temp*temp*dt_M;
            }

            sum_x += org_vx[chosen[i]];
            sum_y += org_vy[chosen[i]];
            if(sum_x < 0 || sum_y < 0){
                return;
            }
        }

    if(sum_x == 0 && sum_y == 0){
        vector<pair<int, int> > comb;
        sum_weight += weight_temp;
        for(int i=0;i<k;i++){
            comb.push_back(make_pair(0, i));
        }
        
        for(int i=0;i<n;i++){
            comb[chosen[i]].first += 1;
        }
        // sort
        // sort(comb.begin(), comb.end(), greater<pair<int, int> >());

        //output
        // for(int i=0;i<n;i++){
        //     cout << chosen[i] << ' ';
        // }
        //     cout << chosen[n-1];
        // cout << ' ';
        // for(int i=0;i<r;i++){
        //     cout << comb[i].second << ' ';
        // }
        // cout << endl;
        // cnt++;
        }
        return;
    }

    for (int i = 0; i < k; i++)
    {
        chosen[index] = i;
        Exhaustive(chosen, index + 1, n, 0, end, weight);
    }
    return;
}

void ExhaustiveRepetition(int k, int n, vector<double>& weight)
{
    int chosen[n+1];
    Exhaustive(chosen, 0, n, 0, k, weight);
}