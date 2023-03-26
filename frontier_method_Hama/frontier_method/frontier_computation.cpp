#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <utility>
#include <set>

using namespace std;

#define TRUE 1
#define FALSE 0
#define MAX_ID 10000000
const int k = 10; //状態遷移数
const int n = 5; //ステップ数
int org_vx[k+1] = {  -2, -1, -1,  0,  0, 0,  1,  1,  2,  2};
int org_vy[k+1] = {  0,  0,  1, -2, -1,  0, -1,  0, -1,  0};
vector<int> vx(org_vx, org_vx+k+1);
vector<int> vy(org_vy, org_vy+k+1);
const int sx = 1;
const int sy = 1;
const int gx = 0;
const int gy = 0;
int cnt = 0;
double weight=0;

void Combination(int chosen[], int index, int n, int start, int end, vector<vector<pair<int, int> > >& comb);
void CombinationRepetition(int k, int n, vector<vector<pair<int, int> > >& comb);
void build_zdd(vector<vector<int> >& zdd, vector<vector<int> >& edge, vector<int>& step_num, vector<int>& step_head);
double weight_computation(vector<vector<int> >& zdd, vector<vector<int> >& edge, vector<int>& step_num, vector<int>& step_head, vector<vector<pair<int, int> > >& comb);
void comp_next(int step, int next_node, int comp, int tempx, int tempy, vector<vector<int> >& zdd, vector<vector<int> >& edge, vector<int>& step_num, vector<int>& step_head, vector<vector<pair<int, int> > >& comb, vector<int>& vx, vector<int>& vy);

int main(int argc, char *argv[]) {

    chrono::system_clock::time_point start, end;

    start = chrono::system_clock::now();
    vector<vector<int> > zdd,edge;
    vector<int> step_num;
    vector<int> step_head;
    vector<vector<pair<int, int> > > comb;
    zdd.resize(MAX_ID);
    edge.resize(MAX_ID);
    step_num.resize(n+2);
    step_head.resize(n+2);

    //計算時間計測_開始地点
    start = chrono::system_clock::now();

    //到達組合せの探索
    CombinationRepetition(k, n, comb);

    //フロンティア法による木構造の構築
    build_zdd(zdd, edge, step_num, step_head);

    //統計量の計算
    weight = weight_computation(zdd, edge, step_num, step_head, comb);

    //統計量の出力
    cout << std::setprecision(16) << weight << endl;

    //計算時間計測_終了地点
    end = chrono::system_clock::now();

    //計算時間の出力
    double time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    printf("time %lf[ms]\n", time);

    return 0;
}

//到達組合せの探索
void Combination(int chosen[], int index, int n, int start, int end, vector<vector<pair<int, int> > >& comb)
{
    if (index == n)
    {
        int sum_x=sx;
        int sum_y=sy;
        for (int i = 0; i < n; i++){
            sum_x += vx[chosen[i]];
            sum_y += vy[chosen[i]];
        }
        if(sum_x == gx && sum_y == gy){
            comb.push_back(vector<pair<int, int> >());
            for(int i=0;i<k;i++){
                comb.at(comb.size()-1).push_back(make_pair(0, i));
            }
            
            for(int i=0;i<n;i++){
                comb.at(comb.size()-1).at(chosen[i]).first += 1;
            }
            //sort
            sort(comb.at(comb.size()-1).begin(), comb.at(comb.size()-1).end(), greater<pair<int, int> >());

            // output
            // 到達組合せの列挙
            // for(int i=0;i<k;i++){
            //     cout << comb.at(comb.size()-1).at(i).first << ' ';
            // }
            // cout << ' ';
            // for(int i=0;i<k;i++){
            //     cout << comb.at(comb.size()-1).at(i).second << ' ';
            // }
            // cout << endl;
            // cnt++;

        }
        return;
    }

    for (int i = start; i <= end; i++)
    {
        chosen[index] = i;
        Combination(chosen, index + 1, n, i, end, comb);
    }
    return;
}

void CombinationRepetition(int k, int n, vector<vector<pair<int, int> > >& comb)
{
    int chosen[n+2];
    Combination(chosen, 0, n, 0, k-1, comb);
}


//木構造の構築
void build_zdd(vector<vector<int> >& zdd, vector<vector<int> >& edge, vector<int>& step_num, vector<int>& step_head)
{
    chrono::system_clock::time_point start, end;
    int tempx,tempy;
    int cnt=0;

    //分割数による分類のファイル読み込み
    string a = to_string(n);
    string file_name = "partition";
    file_name += a;
    file_name += ".txt";
	ifstream ifs;

    ifs.open(file_name,ios::in);
	if (ifs.fail()) {
	   exit(0);
	}
	string str;
    int ret;
    int temp=0,next;
    int flag=TRUE;
    
    // start = chrono::system_clock::now();
    for(int i=0;i<MAX_ID;i++){
        zdd.at(i).resize(k);
        edge.at(i).resize(k);
    }
    // end = chrono::system_clock::now();

	while (getline(ifs,str)) {
	      stringstream ss(str);
          int i=0;
          flag = TRUE;
	      while (!ss.eof()){
              ss >> ret;
              if(i==k){
                  flag = FALSE;
                  break;
              }
              zdd.at(temp).at(i)=ret;
              i++;
          }
          if(flag == TRUE){
              temp++;
          }else{
              for(int j=0;j<k;j++){
                  zdd.at(temp).at(j)=0;
              }
          }
	}

    step_num.at(0)=temp;
    step_head.at(0)=0;

    for(int i=1;i<n+1;i++){ //ステップ数
        step_head.at(i)=step_head.at(i-1)+step_num.at(i-1);

        for(int f=0;f<step_num.at(i-1);f++){ //1つ前のステップからフロンティアを構築する
            temp = step_head.at(i-1)+f;

            for(int j=0;j<k;j++){ //状態遷移数
                if(zdd.at(temp).at(j) >= 1){ //遷移がある場合
                    step_num.at(i) = step_num.at(i)+1; //次のステップのフロンティア数を増やす
                    next = step_head.at(i)+step_num.at(i)-1;
                    flag = TRUE;

                    zdd.at(next)=zdd.at(temp);
                    zdd.at(next).at(j) = zdd.at(next).at(j)-1; //移動した分の遷移を減らす

                    for(int l=0;l<step_num.at(i)-1;l++){ //重複の確認

                        if(zdd.at(step_head.at(i)+l)==zdd.at(next)){ //重複している場合
                            edge.at(temp).at(j)=l+1; //既にあるフロンティアと繋ぐ
                            flag = FALSE;
                            break;
                        }
                    }

                    if(flag == TRUE){ //重複がなければ
                        edge.at(temp).at(j)=step_num.at(i); //追加したフロンティアと繋ぐ
                    }else{
                        step_num.at(i) = step_num.at(i) - 1; //追加したフロンティアを削除
                    }

                }
            }

        }

    }

    // output
    //木構造の出力
    // for(int i=0;i<n+1;i++){
    //     for(int j=0;j<step_num.at(i);j++){
    //         cout << setw(2) << j+1 << ": ";
    //         for(int l=0;l<zdd.at(step_head.at(i)+j).size();l++){
    //             cout << zdd.at(step_head.at(i)+j).at(l) << ' ';
    //         }
    //         cout << ' ';

    //         for(int l=0;l<edge.at(step_head.at(i)+j).size();l++){
    //             cout << edge.at(step_head.at(i)+j).at(l) << ' ';
    //         }
    //         cout << endl;
    //     }
    //     cout << endl;
    //     cnt+=step_num.at(i);
    // }

    return;
}

double weight_computation(vector<vector<int> >& zdd, vector<vector<int> >& edge, vector<int>& step_num, vector<int>& step_head, vector<vector<pair<int, int> > >& comb)
{
    chrono::system_clock::time_point start, end;
    double dt = 0.1;
    double epsilon = 1.0;
    double n_11 = 0.5;
    double n_22 = 0.5;
    double x_1 = 1.0;
    double x_2 = 0.5;
    double sum_all = 0;

    double dt_M = dt / (double)n;

    const int r = 14;
    double weight[r+1];
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

    for(int i=0;i<r;i++){
        for(int j=0;j<k;j++){
            if(walk_x[i]==vx[j] && walk_y[i]==vy[j]){
                walk_corres[i] = j;

            }
        }
    }

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

    //
    vector<int> comb_number;
    comb_number.resize(comb.size());

    for(int i=0;i<comb.size();i++){
        for(int j=0;j<step_num.at(0);j++){
            int flag = TRUE;
            for(int l=0;l<k;l++){
                if(comb.at(i).at(l).first != zdd.at(j).at(l)){
                    flag = FALSE;
                    break;
                }
            }
            if(flag==TRUE){
                comb_number.at(i) = j;
            }
        }
    }

    vector<int> x_coordinate,y_coordinate;
    vector<double> weight_temp;
    int now,next;
    double temp;
    vector<vector<int> > next_array;
    next_array.resize(n+1);
    weight_temp.resize(zdd.size(), 0);
    x_coordinate.resize(zdd.size(), 0);
    y_coordinate.resize(zdd.size(), 0);

    for(int i=0;i<comb.size();i++){ //到達組合せ
        fill(weight_temp.begin(), weight_temp.end(), 0);
        fill(x_coordinate.begin(), x_coordinate.end(), 0);
        fill(y_coordinate.begin(), y_coordinate.end(), 0);
        weight_temp.at(comb_number.at(i))=1;
        x_coordinate.at(comb_number.at(i)) = sx;
        y_coordinate.at(comb_number.at(i)) = sy;
        // cout << comb.at(i).at(0).second << endl; 代入する状態遷移番号
        for(int j=0;j<n+1;j++){ //ステップ数
            for(int l=0;l<step_num.at(j);l++){ //フロンティアのノード
                now = step_head.at(j)+l;
                // cout << now << ": " << step_num.at(j) << ' ' << weight_temp.at(now) << ' ' << x_coordinate.at(now) << ' ' << y_coordinate.at(now) << endl;
                if(x_coordinate.at(now) < 0 || y_coordinate.at(now) < 0){
                    continue;
                }else if(weight_temp.at(now) != 0){
                    for(int m=0;m<k;m++){ //状態遷移数
                        if(edge.at(now).at(m) != 0){
                            next = step_head.at(j+1) + edge.at(now).at(m) - 1;

                            if(comb.at(i).at(m).second == 0){
                                temp = weight[0]*x_coordinate.at(now)*(x_coordinate.at(now)-1);
                            }else if(comb.at(i).at(m).second == 1){
                                temp = weight[1]*x_coordinate.at(now);

                            }else if(comb.at(i).at(m).second == 2){
                                temp = weight[2]*x_coordinate.at(now);

                            }else if(comb.at(i).at(m).second == 3){
                                temp = weight[3]*y_coordinate.at(now)*(y_coordinate.at(now)-1);

                            }else if(comb.at(i).at(m).second == 4){
                                temp = weight[4]*y_coordinate.at(now);
                                temp += weight[5]*y_coordinate.at(now);
                                temp += weight[6]*y_coordinate.at(now);

                            }else if(comb.at(i).at(m).second == 5){
                                temp = weight[7]*y_coordinate.at(now);
                                temp += weight[8]*y_coordinate.at(now);

                            }else if(comb.at(i).at(m).second == 6){
                                temp = weight[9]*y_coordinate.at(now);
                                temp += weight[10]*y_coordinate.at(now);

                            }else if(comb.at(i).at(m).second == 7){
                                temp = weight[11]*y_coordinate.at(now);

                            }else if(comb.at(i).at(m).second == 8){
                                temp = weight[12]*y_coordinate.at(now);

                            }else if(comb.at(i).at(m).second == 9){
                                temp = weight[13]*y_coordinate.at(now);

                            }


                            // for(int n1_i=0;n1_i<index_rate1[m];n1_i++){
                            //     temp = temp*(x_coordinate.at(now) - i);
                            // }
                            // for(int n2_i=0;n2_i<index_rate2[m];n2_i++){
                            //     temp = temp*(y_coordinate.at(now) - i);
                            // }
                            if(comb.at(i).at(m).second == 5){
                                weight_temp.at(next) =  weight_temp.at(next) + weight_temp.at(now)*(1+temp*dt_M);
                            }else{
                                weight_temp.at(next) =  weight_temp.at(next) + weight_temp.at(now)*temp*dt_M;
                            }
                            x_coordinate.at(next) = x_coordinate.at(now) + vx[comb.at(i).at(m).second];
                            y_coordinate.at(next) = y_coordinate.at(now) + vy[comb.at(i).at(m).second];
                        }
                    }
                }
            }
        }
        sum_all += weight_temp.at(next);
    }

    return sum_all;
}