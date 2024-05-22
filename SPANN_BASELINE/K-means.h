//
// Created by 郑荘霖 on 2024/5/7.
//
#include<bits/stdc++.h>
#include <string>
#include "Point.h"

using namespace std;

#ifndef SPANN_BASELINE_K_MEANS_H
#define SPANN_BASELINE_K_MEANS_H
////////////////////////////////
float L = 2; // 聚类平衡惩罚因子  L up  |  time down  |  recall down
float C = 1; // 重分配临界因子
float Q = 1.15; // 查询临界因子
int _k_ = 2; // 重分配数量
int query = 0;  //候选点数量
////////////////////////////////

class K_means {
public:
    int k = 5;
    int N = 10;
    int l = 5;
    vector<Point> centroids;
    vector<vector<Point>> clusters;
    vector<Point> points;
    vector<K_means *> children;
    vector<unordered_set<int>> cluster_ids;
    set<int> res;

    /*
     * 构造函数
     * O(k), k为每层聚类数量
     */
    K_means(vector<Point> &points, int k, int N, int l) {
        this->k = k;
        this->N = N;
        this->l = l;
        this->points = points;
        this->centroids.resize(k);
        this->clusters.resize(k);
        this->children.resize(k);
        this->cluster_ids.resize(k);
    }

    /*
     * L2距离
     * O(m), m为向量维度
     */
    double L2_distance(const Point &v1, const Point &v2) {
        float sum = 0;
        for (int i = 0; i < v1.coordinates.size(); i++) {
            sum += pow(v1.coordinates[i] - v2.coordinates[i], 2);
        }
        return sqrt(sum);
    }

    /*
     * 计算质心
     * O(n*m), n为聚类大小, m为向量维度
     */
    Point get_centroid(const vector<Point> &cluster) {
        if (cluster.empty()) return Point{};
        size_t coordSize = cluster[0].coordinates.size();
        Point centroid;
        centroid.coordinates.resize(coordSize, 0.0f);
        for (size_t i = 0; i < coordSize; ++i) {
            float sum = 0.0f;
            for (const auto &point: cluster) {
                sum += point.coordinates[i];
            }
            // 计算平均值并存储
            centroid.coordinates[i] = sum / cluster.size();
        }
        /**
         * 将质心转化为距离其最近的点
         */
        double min_distance = numeric_limits<double>::max();
        Point min_centroid_point;
        for (int i = 0; i < cluster.size(); ++i) {
            double l2Distance = L2_distance(centroid, cluster[i]);
            if (l2Distance < min_distance) {
                min_distance = l2Distance;
                min_centroid_point = cluster[i];
            }
        }
        centroid = min_centroid_point;
        return centroid;
    }

    /*
     * 计算points所属的cluster聚类
     * O(n*m*k), n是向量个数, m是向量维度, k为聚类数量
     */
    void assign_clusters() {
        for (const auto &point: this->points) {
            double min_distance = numeric_limits<double>::max();  // 初始化最小距离为正无穷
            int cluster_index = 0;
            vector<double> point_distances;
            for (int i = 0; i < k; ++i) {
                /**
                 * 1.找到最近的质心
                 * 2.加入修正项，使得每个聚类的点个数尽量平均
                 * 3.L为修正因子
                 */
                double dist = L2_distance(point, this->centroids[i]);
                point_distances.push_back(dist);
                double penalty = L * this->clusters[i].size();
                dist += penalty;
                if (dist < min_distance) {
                    min_distance = dist;
                    cluster_index = i;
                }
            }
            this->clusters[cluster_index].push_back(point);
            this->cluster_ids[cluster_index].insert(point.point_id);
            duplicate_points(point_distances, cluster_index, point_distances[cluster_index], point);
        }
    }

    /*
     * 多次分配同一向量
     * O(k*m), k为聚类数量, m是向量维度
     */
    void duplicate_points(vector<double> &point_distances, int &cluster_index, double &min_distance, Point point) {
        int duplicate = 0;
        vector<int> duplicate_index;
        duplicate_index.push_back(cluster_index);
        int random = rand();
        for (int i = 0; i < point_distances.size(); ++i) {
            int item = (i + random) % point_distances.size();
            bool skip = false;
            bool copy = false;
            /**
             * 如果当前点与另一个质心距离小于min_distance * C，则将当前点加入另一个质心所属的聚类中
             * 能够较好提高召回率
             * min_distance应该为原本的dist
             */
            if (item != cluster_index && point_distances[item] < min_distance * C) {
                copy = true;
            }
            if(!copy) {
                continue;
            }
            for (int j = 0; j < duplicate_index.size(); ++j) {
                if (L2_distance(centroids[duplicate_index[j]], centroids[item]) < point_distances[item]) {
                    skip = true;
                    break;
                }
            }
            if (!skip) {
                duplicate_index.push_back(item);
                if (this->cluster_ids[item].find(point.point_id) == this->cluster_ids[item].end()) {
                    this->clusters[item].push_back(point);
                }
                this->cluster_ids[item].insert(point.point_id);
                duplicate++;
                if (duplicate >= _k_) break;
            }
        }
    }

    /*
     * 初始化聚类点起
     * O(k),k为聚类数量
     */
    void initializeCentroids() {
        for (int i = 0; i < this->k; ++i) {
            int randomIndex = rand() % this->points.size();
            this->centroids[i] = this->points[randomIndex];
        }
    }


    // 运行K-means
    void run(int max_iterations) {
        initializeCentroids();
        // 初始化clusters O(n*m*k)
        this->assign_clusters();
        /*
         * 开始迭代
         * O(i*n*m*k) i为迭代次数
         */
        for (int iter = 0; iter < max_iterations; ++iter) {
            // 更新centroids
            vector<Point> new_centroids(this->k);
            for (int i = 0; i < this->k; ++i) {
                new_centroids[i] = get_centroid(this->clusters[i]);
            }
            // 清除之前记录的聚类信息
            for (auto &cluster: this->clusters) {
                cluster.clear();
            }
            for (auto &ids: this->cluster_ids) {
                ids.clear();
            }
            // 重新计算质点
            this->assign_clusters();
            // 判断是否终止迭代
            bool converged = true;
            for (int i = 0; i < this->k; ++i) {
                if (this->centroids[i].point_id != new_centroids[i].point_id) {
                    converged = false;
                    break;
                }
            }
            if (converged) {
                break;
            }
            this->centroids = new_centroids;
        }
        this->check_clusters();
    }

    /*
     * 检验clusters的大小, 如果大于l,则进行子聚类
     * O(n*m*k*log_k(N))
     */
    void check_clusters() {
        for (int i = 0; i < this->clusters.size(); ++i) {
            vector<Point> cluster = this->clusters[i];
            if (cluster.size() > this->l) {
                /**
                 * 修正k的大小
                 */
                int curr_k = cluster.size() / this->l + 1;
                if (this->k < curr_k) {
                    curr_k = this->k;
                }
                K_means *child_cluster = new K_means(cluster, curr_k, this->N, this->l);
                child_cluster->run(10);
                this->children[i] = child_cluster;
            } else {
                this->children[i] = nullptr;
            }
        }
    }

    /*
     * 输出clusters中各个cluster的大小
     * O(k)
     */
    void print_clusters() {
        for (int i = 0; i < this->k; ++i) {
            cout << "Cluster " << i << ": " << this->clusters[i].size() << endl;
        }
    }

    // 查询目标向量的ANN
    void query_ANN(const Point &query_point, int num) {
        priority_queue<pair<float, int>, vector<pair<float, int>>, less<>> candidates;
        unordered_set<int> visited;
        candidate_points(query_point, num, candidates, visited);
        res.clear();
        while (!candidates.empty()) {
            res.insert(candidates.top().second);
            candidates.pop();
        }
    }

    /*
     * 获取候选点集
     */
    void candidate_points(const Point &query_point, int num,
                          priority_queue<pair<float, int>, vector<pair<float, int>>, less<>> &candidates,
                          unordered_set<int> &visited) {
        // 遍历获取距离
        priority_queue<pair<float, int>, vector<pair<float, int>>, less<>> distances;
        double min_distance = numeric_limits<double>::max();
        for (size_t i = 0; i < this->centroids.size(); i++) { //O(k*m)
            float distance = this->L2_distance(query_point, this->centroids[i]);
            if (distance < min_distance) {
                min_distance = distance;
            }
            // 利用priority_queue存储距离最小的num个向量
            if (distances.size() < num) {
                distances.push(make_pair(distance, i));
            } else if (distance < distances.top().first) {
                distances.pop();
                distances.push(make_pair(distance, i));
            }
        }
        while (!distances.empty()) {
            /**
             * 如果当前点与另一个质心距离小于min_distance*Q，则将当前点加入候选点集
             */
            vector<Point> &vectors = this->clusters[distances.top().second];
            if (distances.top().first < min_distance * Q) {
                if (this->children[distances.top().second] != nullptr) {
                    this->children[distances.top().second]->candidate_points(query_point, num, candidates, visited);
                } else {
                    for (auto &point: vectors) {
                        if (visited.find(point.point_id) == visited.end()) {
                            query++;
                            if (candidates.size() < num) {
                                candidates.push(make_pair(this->L2_distance(query_point, point), point.point_id));
                            } else if (this->L2_distance(query_point, point) < candidates.top().first) {
                                candidates.pop();
                                candidates.push(make_pair(this->L2_distance(query_point, point), point.point_id));
                            }
                            visited.insert(point.point_id);
                        }
                    }
                }
            }
            distances.pop();
        }
    }

    /*
     * 输出vector向量
     * O(m)
     */
    vector<float> print_vector(vector<float> v, bool out = true) {
        sort(v.begin(), v.end());
        if (out) {
            for (auto i: v) {
                cout << i << " ";
            }
            cout << "\n";
        }
        return v;
    }

    void print_size() {
        // 输出clusters中各个cluster的大小
        for (int i = 0; i < this->k; ++i) {
            cout << "Cluster " << i << ": " << this->clusters[i].size() << endl;
        }
    }
};


#endif //SPANN_BASELINE_K_MEANS_H
