//
// Created by 郑荘霖 on 2024/5/10.
//
#include<bits/stdc++.h>
#include <string>

using namespace std;

#ifndef SPANN_BASELINE_POINT_H
#define SPANN_BASELINE_POINT_H
struct Point {
    vector<float> coordinates;
    int point_id = -1;
};
void readDataFromFile(const string &filename, vector<Point> &points, int dataset);
void txt_read(const string &filename, vector<Point> &points, int dataset);
void ivecs_read(const string &fname, vector<Point> &points, int dataset);

void readDataFromFile(const string &filename, vector<Point> &points, int dataset){
    vector<string> parts;
    size_t start = 0, end;
    while ((end = filename.find('.', start)) != string::npos) {
        parts.push_back(filename.substr(start, end - start));
        start = end + 1;
    }
    if (start < filename.length()) {
        parts.push_back(filename.substr(start));
    }
    //获取parts的最后一个数据
    string suffix = parts.back();
    if (suffix == "txt") {
        txt_read(filename, points, dataset);
    } else if (suffix == "ivecs" || suffix == "fvecs") {
        ivecs_read(filename, points, dataset);
    } else {
        cerr << "不支持的文件格式: " << filename << std::endl;
    }
}

void txt_read(const string &filename, vector<Point> &points, int dataset) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "无法打开文件: " << filename << std::endl;
        return;
    }
    string line;
    int index = 0;
    while (getline(file, line)) {
        Point point;
        istringstream iss(line);
        float value;
        int length = 0;
        while (iss >> value) {
            point.coordinates.push_back(value);
            length++;
            if (length == 128) break;
        }
        point.point_id = index;
        points.push_back(point);
        index++;
        if (index == dataset) break;
    }
    file.close();
}

void ivecs_read(const string &fname, vector<Point> &points, int dataset) {
    const size_t VECTOR_DIM = 128;
    ifstream file(fname, ios::binary);
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << fname << std::endl;
        return;
    }
    int index = 0;
    vector<float> vector(VECTOR_DIM);
    while (file.read(reinterpret_cast<char *>(vector.data()), VECTOR_DIM * sizeof(float))) {
        Point point;
        point.coordinates = vector;
        point.point_id = index;
        points.push_back(point);
        index++;
        if (index == dataset) break;
    }
    file.close();
}


#endif //SPANN_BASELINE_POINT_H
