#include <iostream>
#include <vector>
#include <tuple>
//#include "IO.h"
//#include "hypercube.h"
#include "subann.hpp"
#include "config.hpp"
#include "curve.hpp"
#include "point.hpp"
#include "frechet.hpp"
#include <unordered_map>
#include "utils.hpp"
#include "timeseries.hpp"
#include <ctime>
#include <ratio>
#include <chrono>
#include <cstdlib>
#include <ctime>
#define N 10000

namespace fc = Frechet::Continuous;
using namespace std;

int main() {
	int i;
    std::srand(std::time(nullptr));


    vector<double> signal = {1.0, 2.0, 3.0, 2.5, 4.0, 5.0, 4.5, 3.0, 2.0};
	vector<vector<coordinate_t>> signal2d = {{1.0}, {2.0}, {3.0}, {2.5}, {4.0}, {5.0}, {4.5}, {3.0}, {2.0}};
    vector<vector<coordinate_t>> signal2d2 = {{1.0}, {2.0}, {3.0}, {2.5}, {4.0}, {5.0}, {4.5}, {3.0}, {-1.0}};

    //test signature
	double R = 0.1;
    vector<int> indx = signature(signal, R);
    for (i = 0; i < indx.size(); i++) {
        cout << indx[i] << " ";
    }
	cout<<endl;
	//test canonical
	vector<int> indx2 = canonical(signal);
    for (i = 0; i < indx2.size(); i++) {
        cout << indx2[i] << " ";
    }
	cout<<endl;

	//test remove reduntant
	vector<double> new_signal = remove_redundant(signal);
    for (int i = 0; i < new_signal.size(); i++) {
        cout << new_signal[i] << " ";
    }
    cout<<endl;

	//test get_sequence_id
	double w = 0.5;
    double s = static_cast<double>(std::rand())*w / RAND_MAX;

	vector<int> frechet_id = get_frechet_id(signal, w, s);
    cout << "Sequence ID: ";
    for (int i = 0; i < frechet_id.size(); i++) {
        cout << frechet_id[i] << " ";
    }
    cout << endl;

	//test weak_frechet_id
	//vector<tuple<int, int>> result = weak_frechet_id(signal, w, s);
    vector<int> weak_frechet_id = get_weak_frechet_id(signal, w, s);
    for (int i = 0; i < weak_frechet_id.size(); i++) {
        //cout << "(" << get<0>(result[i]) << ", " << get<1>(result[i]) << ") ";
        cout << weak_frechet_id[i];
    }
    cout << endl;

	//test frechet computation
    //Curve  signal_curve1(signal2d);
    // Curve  signal_curve2(signal2d2);

    // cout<<"First curve:"<<signal_curve1<<endl;
    // cout<<"Second curve:"<<signal_curve2<<endl;

    // cout<<"Distance between two curves:"<<fc::distance(signal_curve1, signal_curve2).value<<endl;

	//concatenate weak frechet id and frechet id
    // std::vector<int> id;
    // id = get_id(signal, w, s); 
    // for (int i = 0; i < id.size(); i++) {
    //     cout << id[i];
    // }
    // cout << endl;

    //unordered_map<vector<int>, vector<pair<int,int>>, VectorHasher> map;
    // // Insert an element into the map using the vector as the key
    // vector<int> list{1,2,3,42};
    // map.emplace(id, list);

    // // Check if the vector is a key in the map
    // if (map.count(id)) {
    //     std::cout << "id is a key in the map" << std::endl;
    // }
    // std::cout<<map[id][3]<<std::endl;
    
    // //Preprocessing
    // vector<pair<int, int>> combinations;
    // // Create all combinations with replacement of the signal indices
    // for (int i = 0; i < signal.size(); i++) {
    //     for (int j = i; j < signal.size(); j++) {
    //         combinations.push_back(make_pair(i, j));
    //     }
    // }
    // // Loop through each combination and store corresponding subsequence
    // for (auto combination : combinations) {
    //     //cout << "(" << get<0>(combination) << ", " << get<1>(combination) << ") ";
    //     auto start_it = signal.begin() + get<0>(combination);
    //     auto end_it = signal.begin() + get<1>(combination)+1;
    //     // Create new vector with the desired subsequence
    //     vector<double> sub_signal(start_it, end_it);
        
    //     id = get_id(sub_signal, w, s);
    //     cout<<endl;
    //     if (map.count(id)){
    //         map[id].push_back(combination);
    //     }
    //     else{
    //         vector<pair<int,int>> newlist {combination};
    //         map.emplace(id, newlist);
    //     }

    // }
    
    //TODO: frechet.cpp: add cerr, 


    //vector<vector<coordinate_t>> signal1 = {{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
    //vector<vector<coordinate_t>> signal2 = {{1.0, 1.0}, {2.0, 4.0}, {5.0, 6.0}};

    vector<vector<double>> signal1 = {{1.0}, {3.0}, {6.0}, {1.0}, {3.0}, {6.0}, {1.0}, {3.0}, {6.0}, {1.0}, {3.0}, {6.0},{1.0}, {3.0}, {6.0}, {1.0}, {3.0}, {6.0}, {1.0}, {3.0}, {6.0}, {1.0}, {3.0}, {6.0},{1.0}, {3.0}, {6.0}, {1.0}, {3.0}, {6.0}, {1.0}, {3.0}, {6.0}, {1.0}, {3.0}, {6.0},{1.0}, {3.0}, {6.0}, {1.0}, {3.0}, {6.0}, {1.0}, {3.0}, {6.0}, {1.0}, {3.0}, {6.0},{1.0}, {3.0}, {6.0}, {1.0}, {3.0}, {6.0}, {1.0}, {3.0}, {6.0}, {1.0}, {3.0}, {6.0}};
    // for(int i = 0;i < 100;i++){
    //     signal1.push_back(signal1[1]);
    // }
    vector<vector<coordinate_t>> signal2 = {{1.0}, {4.0}, {6.0}};

	TimeSeries ts(signal1);
    TimeSeries query(signal2);
    // cout<<ts.size()<<endl;

    for(auto a: ts.get_index()){
        cout<<a<<" ";
    }        
    cout<<endl;
    //SubAnnDS ds(ts, 10, 5,2);

    SubAnnDS ds(ts, 10, 5);


    // cout<<ds<<endl;
    // vector<int> key{-1,-1};
    cout<<"AAAAAAAAAAA"<<endl;
    //std::vector<std::pair<int,int>>* answers = ds.range_query(query, 100);
    std::vector<std::pair<int,int>> answers = ds.range_query(query, 0.0001);
    //std::vector<int> answers = ds.range_query(query, 100000);
    cout<<answers.size()<<endl;
    ds.stats();
    // cout<<"Reporting answers"<<endl;
    // for(auto x: answers){
    //     cout<<"("<<x.first<<", "<<x.second<<")"<<endl;
    // }
    // cout<<ds;






    //TODOs: 
    //ann query: what happens if bucket is empty8
    //type of index in TimeSeries
    // offline metching
    //paralellize constructor of SubANN
    //range_query: return curves? 
    // fix use of either timeseries or curve
    // figure out what's going on with ostream& operator<<(ostream& os, const TimeSeries& obj);
    // weights w's
    // parallelization: gomp
    // wrapper to signature 
	// templates
    
   
}






