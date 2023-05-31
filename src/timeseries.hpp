#ifndef TIMESERIES_H
#define TIMESERIES_H
#include <iostream>
#include <vector>

using namespace std;

class TimeSeries {
private:
    vector<vector<double>> data;
    vector<int> index;
    int dim;
public:
    // constructor
    inline TimeSeries(){};
    TimeSeries(const vector<vector<double>>& data);
    TimeSeries(const vector<vector<double>>& data, const vector<int>& index);
    
    // get a specific point in the time series
    vector<double> operator[](int i) const {
        return data[i];
    }
    vector<vector<double>> values() const {
        return data;
    }
    vector<int> get_index() const{
        return index;
    }
    const TimeSeries subseries(int i, int j) {
        if((i>j)||(i<0)||(j>=data.size())){
            cerr<<"Subseries: wrong indices"<<endl;
        }
        //cout<<"subseries "<<i<<"  "<<j<<endl;
        //cout<<data.size()<<endl;
        std::vector<std::vector<double>> sub_v(data.begin()+i, data.begin()+j+1);
        //cout<<"subseries "<<endl;
        vector<int> sub_index(index.begin()+i, index.begin()+j+1);
        // for (std::vector<std::vector<double>>::iterator it = begin; it != end; ++it) {
        //     sub_v.push_back(*it);
        // }
        return TimeSeries(sub_v, sub_index);
    }
    TimeSeries iloc(vector<int> indices){
        vector<vector<double>> subseq;
        for(int i: indices){
            subseq.push_back(data[i]);
        }
        return TimeSeries(subseq);
    }
    const int get_dim() const {
        return dim;
    }
    // get the number of points in the time series
    int size() const {
        return data.size();
    }
    std::vector<vector<double>>::const_iterator begin()const {
        return data.begin();
    }
    std::vector<vector<double>>::const_iterator end()const {
        return data.end();
    }
    //project to a certain axis
    vector<double>  projectAxis(int j);
    friend ostream& operator<<(ostream& os, const TimeSeries& obj);
};
#endif
