#include "K-means.h"
#include "Baseline.h"

#ifdef __linux__
#include "monitor.h"
#endif

#ifdef __linux__
program_t begProg, endProg;
#endif

// int argc, char *argv[]
int main(int argc, char *argv[]) {
    // 读入参数
    string parameter_file = argv[1];
    string prefix = argv[2];
    ifstream file_arg(parameter_file);
    if (!file_arg.is_open()) {
        std::cerr << "Error opening file: " << parameter_file << std::endl;
        return 1;
    }
    vector<float> parameters;
    float parameter;
    while (file_arg >> parameter) {
        parameters.push_back(parameter);
    }
    int k = int(parameters[0]);
    int l = int(parameters[1]);
    int num = int(parameters[2]);
    int dataset = int(parameters[3]);
    int times = int(parameters[4]);
    L = int(parameters[5]);
    C = parameters[6];
    Q = parameters[7];
    ///////////////////////////////////
    vector<Point> points;
    vector<Point> targets;
    vector<set<int>> SPANN_res, baseline_res;
    string filename = prefix + "/data/SIFT1M/sift_base.fvecs";
    string target = prefix + "/data/SIFT1M/sift_query.fvecs";
    /*
     * k 层聚类数量
     * N 总聚类数量
     * l 每个聚类的上限数量
     * num ANN中最临近向量数量
     * dataset 数据库大小
     */
    int N = 100;
    readDataFromFile(filename, points, dataset);
    readDataFromFile(target, targets, times);
    // k-means
    K_means k_means(points, k, N, l);
    k_means.run(10);
#ifdef __linux__
    save_time(begProg);
#endif
    int index = 0;
    while (index < times) {
        k_means.query_ANN(targets[index], num);
        set<int> copy = set<int>(k_means.res);
        SPANN_res.push_back(copy);
        index++;
    }
#ifdef __linux__
    save_time(endProg);
    double usedTime = calc_time(begProg, endProg);
    printf("SPANN: %.6lf \n", usedTime);
#endif
    // baseline
    Baseline baseline(points);
#ifdef __linux__
    save_time(begProg);
#endif
    index = 0;
    while (index < times) {
        baseline.get_nearest_neighbors_priority_queue(targets[index], num);
        set<int> copy = set<int>(baseline.res);
        baseline_res.push_back(copy);
        index++;
    }
#ifdef __linux__
    save_time(endProg);
    usedTime = calc_time(begProg, endProg);
    printf("priority query: %.6lf \n", usedTime);
#endif
    float sum = 0;
    float sum2 = 0;
    for (int i = 0; i < times; ++i) {
        set<int> s1 = SPANN_res[i];
        set<int> s2 = baseline_res[i];
        set<int> intersection;
        set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                              inserter(intersection, intersection.begin()));
        sum += intersection.size();
        sum2 += s1.size();
    }
    cout << "recall:" << sum / sum2 << endl;
    return 0;
}

/*
string parameter_file = argv[1];
string prefix = argv[2];
ifstream file_arg(parameter_file);
if (!file_arg.is_open()) {
    std::cerr << "Error opening file: " << parameter_file << std::endl;
    return 1;
}
vector<int> parameters;
int parameter;
while (file_arg >> parameter) {
    parameters.push_back(parameter);
}
int k = parameters[0];
int l = parameters[1];
int num = parameters[2];
int dataset = parameters[3];
int times = parameters[4];
*/