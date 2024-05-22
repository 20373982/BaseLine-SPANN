//
// Created by 郑荘霖 on 2024/5/10.
//
#include<bits/stdc++.h>
#include <string>
#include "Point.h"

using namespace std;
#ifndef SPANN_BASELINE_BASELINE_H
#define SPANN_BASELINE_BASELINE_H

class Baseline {
public:
    vector<Point> points;
    set<int> res;

    Baseline(vector<Point> &points) {
        this->points = points;
    }

    double L2_distance(const Point &v1, const Point &v2) {
        float sum = 0;
        for (int i = 0; i < v1.coordinates.size(); i++) {
            sum += pow(v1.coordinates[i] - v2.coordinates[i], 2);
        }
        return sqrt(sum);
    }

    void get_nearest_neighbors_priority_queue(const Point &query_point, int num) {
        // 遍历vectors获取距离
        priority_queue<pair<float, int>, vector<pair<float, int>>, less<>> distances;
        for (size_t i = 0; i < this->points.size(); i++) {
            float distance = L2_distance(query_point, this->points[i]);
            // 利用priority_queue存储距离最小的k个向量
            if (distances.size() < num) {
                distances.push(make_pair(distance, i));
            } else if (distance < distances.top().first) {
                distances.pop();
                distances.push(make_pair(distance, i));
            }
        }
        // 遍历distances输出索引对应的向量
        res.clear();
        while (!distances.empty()) {
            res.insert(distances.top().second);
            distances.pop();
        }
    }

    // 输出vector向量
    vector<float> print_vector(vector<float> v,bool out = true) {
        sort(v.begin(), v.end());
        if(out) {
            for (auto i: v) {
                cout << i << " ";
            }
            cout << "\n";
        }
        return v;
    }
};

#endif //SPANN_BASELINE_BASELINE_H
