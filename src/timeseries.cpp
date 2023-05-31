#include "timeseries.hpp"

/*
* Constructor of the class TimeSeries. 
* @param x The input time series in form of a vector of vectors of doubles.
*/
TimeSeries::TimeSeries(const vector<vector<double>>& x){
    vector<int> temp(x.size());
    for(int i=0; i<x.size();i++){
        temp[i] = i;
    }
    index = temp;
    
    data = x;
    dim = data[0].size();
    for(auto y: data){
        if(y.size()!= dim){
            cerr<<"All points in time series must be of the same dimension."<<endl;
        }
    }
}

/*
* Constructor of the class TimeSeries. 
* @param x The input time series in form of a vector of vectors of doubles.
* @param y The index of the time series in the form of a vector of ints.  
*/
TimeSeries::TimeSeries(const vector<vector<double>>& x, const vector<int>& y){
    if(x.size()!= y.size()){
        cerr<<"Timeseries index does not coincide with its length"<<endl;
    }
    for(auto z: x){
        if(z.size()!= x[0].size()){
            cerr<<"All points in time series must be of the same dimension."<<endl;
        }
    }

    if(y.size()>0){
        for(int i=0; i<(y.size()-1);i++){
            if(y[i]>y[i+1]){
                cerr<<"Index must be monotone"<<endl;
            }
        }
    }    
    data = x;
    index = y;
    dim = data[0].size();
    
    
}

/*
* Projects a multivariate time series to a given axis. 
* @param j Defines the axis to project.
* @return The projected univariate time series in the form of a vector of doubles
*/
vector<double> TimeSeries::projectAxis(int j) {
    vector<double> projection;
    for (int i = 0; i < data.size(); i++) {
        projection.push_back(data[i][j]);  
    }
    return projection;
}    

/*
*Outputs the timeseries to the 'os' stream
*/
ostream& operator<<(ostream& os, const TimeSeries obj) {
    int dim = obj.get_dim();
    for (int i = 0; i < obj.size(); i++) {
        os << "(";
        for (int j = 0; j < dim; j++){
            os << obj[i][j];
            if (j<dim-1){
                os<<", ";
            }
        }
        os << ")";
        if (i<obj.size()-1){
            os<<", ";
        }
    }    
    os <<endl;
    return os;
}