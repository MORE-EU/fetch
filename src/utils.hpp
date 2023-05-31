#ifndef UTILS_H
#define UTILS_H



#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <tuple>
#include <cmath>
#include "timeseries.hpp"

using namespace std;

inline vector<int> signature(vector<double> signal, double R) {
    vector<int> indx = {0};
    bool label = false;
    std::vector<double>::size_type j;
    for (j = 0; j < signal.size(); j++) {
        if (abs(signal[j] - signal[0]) > R) {
            break;
        }
    }
    int b = j;
    for (int i = j; i < signal.size(); i++) {
        if ((signal[i] > max(signal[indx.back()], signal[b]) && signbit(signal[b] - signal[indx.back()]) == 0) || (signal[i] < min(signal[indx.back()], signal[b]) && signbit(signal[b] - signal[indx.back()]) == 1)) {
            b = i;
        }
        else if (abs(signal[i] - signal[b]) > 2*R) {
            indx.push_back(b);
            b = i;
        }
    }
    if (indx.back() < signal.size()-1) {
        indx.push_back(signal.size()-1);
    }
    return indx;
}


inline vector<int> canonical(vector<double> signal) {
    if (signal.size() == 1) {
        return {0};
    }
    int i = 0;
    vector<int> indx = {0};
    //cout<<"size= "<<signal.size()<<" canonical3 "<< ( i < signal.size()-3)<<endl;
    if (signal.size()>=3){
        while ( i < signal.size()-3) {
            int j = i+1;
            int mini = j;
            int maxi = j;
            while (j < signal.size() && j > i && (signal[i] > signal[maxi] || signal[i] < signal[mini])) {
                if (signal[j] > signal[maxi]) {
                    maxi = j;
                }
                if (signal[j] < signal[mini]) {
                    mini = j;
                }
                j = j+1;
            }
            if (abs(signal[i]-signal[maxi]) > abs(signal[i]-signal[mini])) {
                //cout<<"size= "<<signal.size()<<"canonical4 "<<i<<endl;
                indx.push_back(maxi);
                i = maxi;
            }
            else {
                //cout<<"size= "<<signal.size()<<" canonical5 "<<i<<endl;
                indx.push_back(mini);
                i = mini;
            }
            //cout<<"size= "<<signal.size()<<" canonical3 "<<i<<endl;

        }
    }    
    if (indx.back()!= signal.size()-1){
        indx.push_back(signal.size()-1);
    }
    return indx;
}


inline vector<double> remove_redundant(vector<double> signal) {
    // Find the indices of the peaks and valleys in the signal
    vector<int> peaks;
    vector<int> valleys;
    for (int i = 1; i < signal.size()-1; i++) {
        if (signal[i] > signal[i-1] && signal[i] > signal[i+1]) {
            peaks.push_back(i);
        }
        else if (signal[i] < signal[i-1] && signal[i] < signal[i+1]) {
            valleys.push_back(i);
        }
    }
    // Sort the indices and add the first and last elements
    vector<int> inds = peaks;
    inds.insert(inds.begin(), 0);
    inds.push_back(signal.size()-1);
    inds.insert(inds.end(), valleys.begin(), valleys.end());
    sort(inds.begin(), inds.end());
    // Extract the signal values corresponding to the sorted indices
    vector<double> new_signal;
    for (int i = 0; i < inds.size(); i++) {
        new_signal.push_back(signal[inds[i]]);
    }
    return new_signal;
}


inline int get_point_id(double x, double w, double s) {
    return floor((x-s)/w);
}

inline vector<int> get_frechet_id(vector<double> x, double w, double s) {
    vector<int> L;
    for (int i = 0; i < x.size(); i++) {
        if (L.size() == 0 || get_point_id(x[i],w,s) != L.back()) {
            // cout<<"printing in get frechet id"<<endl;
            // cout<<"x: "<<x[i]<<" w: "<<w<<" s: "<<s<<endl;
            L.push_back(get_point_id(x[i],w,s));
        }
    }
    return L;
}
inline vector<int> get_query_frechet_id(vector<double> x, double w, double s, double delta) {
    vector<double> y;
    for(int j: signature(x, delta)){////////////////////// delta MUST BE > r
        y.push_back(x[j]);
    }
    
    std::vector<int> key = get_frechet_id(y, w,  s);
    return key;
}


inline vector< int> get_weak_frechet_id(vector<double>& x, double w, double s) {
    vector<int> L;
    for (auto& p : x) {
        L.push_back(get_point_id(p, w, s));
    }
    vector<int> H;
    for (size_t i = 0; i < L.size() - 1; ++i) {
        
        if (H.empty() ||((L[i] != H[H.size()-2])&& (L[i+1]!=H.back())) &&((L[i] != H.back())&& (L[i+1]!=H[H.size()-2]))){
            H.push_back(L[i]);
            H.push_back(L[i+1]);
        }
    }
    return H;
}

inline vector< int> get_id(TimeSeries series, double w1, double s1, double delta, bool bsign){
    vector<int> frechet_id, weak_frechet_id, id;
    vector<double> x;
    //cout<<"DEBUGGING "<<series.get_dim()<<endl;
    for (int i=0; i<series.get_dim();++i){
        x = series.projectAxis(i);
        vector<double> y;
        if(bsign){
            for(int j: signature(x, delta)){////////////////////// change here THIS MUST BE > r
                y.push_back(x[j]);
            }
        } 
        else{
            y = x;
        }   
        // cout<<endl;
        frechet_id = get_frechet_id(y, w1, s1);
        
        
        id.insert(id.end(), frechet_id.begin(), frechet_id.end());
        //id.insert(id.end(), weak_frechet_id.begin(), weak_frechet_id.end()); // <<<<<<<<<WEAK FRECHET KEY   

    }
    return id;
}

inline vector< int> get_id(TimeSeries series, double w1, double s1, double w2, double s2, double delta, bool bsign){
    /*
    applies weak frechet hashing
    */
    vector<int> frechet_id, weak_frechet_id, id;
    vector<double> x;
    //cout<<"DEBUGGING "<<series.get_dim()<<endl;
    for (int i=0; i<series.get_dim();++i){
        x = series.projectAxis(i);
        vector<double> y;
        if(bsign){
            for(int j: signature(x, delta)){////////////////////// change here THIS MUST BE > r
                y.push_back(x[j]);
            }
        } 
        else{
            y = x;
        }   
        frechet_id = get_frechet_id(y, w1, s1);
        
        vector<double> z;
        for(int j: canonical(x)){////////////////////// change here
            z.push_back(x[j]);
        }
        weak_frechet_id = get_weak_frechet_id(z, w2, s2);


        id.insert(id.end(), frechet_id.begin(), frechet_id.end());
        id.insert(id.end(), weak_frechet_id.begin(), weak_frechet_id.end()); // <<<<<<<<<WEAK FRECHET KEY   

    }
    return id;
}

// Hash function for std::vector<int>
struct VectorHasher {
    std::size_t operator()(const std::vector<int>& v) const {
        std::size_t seed = v.size();
        for (const auto& i : v) {
            seed ^= std::hash<int>()(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

//basic math
inline double dot_product(const vector<double>& v, const vector<double>& u) {
    // returns the dot product of two vectors of the same dimension
    if (v.size()!=u.size()){
        cerr<<"Dot product: Vectors of non-equal dimension"<<endl;
        return std::nan("");
    }
    double result = 0.0;
    for (int i = 0; i < v.size(); i++) {
        result += v[i]*u[i];
    }
    return result;
}

inline double norm(const vector<double>& v) {
    // returns the Euclidean norm of a vector
    return sqrt(dot_product(v,v));
}
inline const std::vector<double> vec_diff(const vector<double>& v, const vector<double>& u){
    if (v.size()!=u.size()){
        cerr<<"Vector difference: Vectors of non-equal dimension"<<endl;
        return {std::nan("")};
    }
    int dim = v.size();
    std::vector<double> diff(dim);
    for (int i = 0; i < dim; i++) {
        diff[i] = v[i] - u[i];
    }
    return diff;
}


inline bool lexcomparePairs(const std::pair<int, int>& a, const std::pair<int, int>& b) {
    // Compare the first elements of the pairs
    if (a.first == b.first){
        return a.second < b.second; 
    }
    return a.first < b.first;
}

#endif