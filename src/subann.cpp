#include "subann.hpp"
#include "frechet.hpp"
#include "curve.hpp"
#include <numeric>
#include <typeinfo>
#include <algorithm>
#include <vector>
#include <typeinfo>
#include <map>
/*
* Constructor of the class SubAnnDS. It implements preprocessing of the data structure
* @param data The input time series.
* @param z1 The parameter defining the side-length of the randomly shifted grid
* @param z2 The parameter defining delta in the delta-signature computations during query
*/
SubAnnDS::SubAnnDS(TimeSeries& data, const double z1, const double z2){
    w1 = z1;
    signal = data;
    delta = z2;
    window = -1;
    s1 = static_cast<double>(std::rand())*w1 / RAND_MAX;
    
    vector<pair<int, int>> combinations;
    // Create all combinations with replacement of the signal indices
    for (int i = 0; i < signal.size(); i++) {
        for (int j = i; j < signal.size(); j++) {
            combinations.push_back(make_pair(i, j));
        }
    }
    // Loop through each combination and store corresponding subsequence
    for (auto combination : combinations) {
        auto start_it = signal.begin() + get<0>(combination);
        auto end_it = signal.begin() + get<1>(combination)+1;
        // Create new vector with the desired subsequence
        vector<vector<double>> sub_vec(start_it, end_it);
        TimeSeries sub_signal(sub_vec);
        vector<int> id = get_id(sub_signal, w1, s1, delta, false);
        if (dict.count(id)){
            dict[id].push_back(combination);
        }
        else{
            vector<pair<int,int>> newlist {combination};
            dict.emplace(id, newlist);
        }
    }
}
/*
* Constructor of the class SubAnnDS. It implements preprocessing of the data structure
* @param data The input time series.
* @param z1 The parameter defining the side-length of the randomly shifted grid
* @param z2 The parameter defining delta in the delta-signature computations during query
* @param window The parameter defining the length of the candidate solutions in the index space
* @param leq Boolean parameter defining if all candidate solutions with length *at most* window will be considered 
*/
SubAnnDS::SubAnnDS(TimeSeries& data, const double z1, const double z2, const double m, const bool lq){
    w1 = z1;
    signal = data;
    delta = z2;
    leq = lq;
    s1 = static_cast<double>(std::rand())*w1 / RAND_MAX;
    int j;    
    window = m;
    std::vector<int> ind =signal.get_index();
    bool bsign = false;
    if(leq){
        for(int i=0; i< data.size();i++){
            j = i;
            while( (j <data.size())&&(ind[j] <= ind[i]+window)){
                vector<vector<double>> sub_vec(signal.begin()+i, signal.begin()+j+1);
                TimeSeries sub_signal(sub_vec);
                vector<int> id = get_id(sub_signal, w1, s1, delta, bsign);
                if (dict.count(id)){
                    dict[id].push_back(std::make_pair(i, j));
                }
                else{
                    vector<pair<int,int>> newlist {std::make_pair(i, j)};
                    dict.emplace(id, newlist);
                }
                j++;
            }
        }
            
    }
    else{
        int n =0;
        for(int i=0; i< data.size();i++){
            auto w_it = std::lower_bound(ind.begin(), ind.end(), ind[i]+window, std::less_equal<int>());
            j = std::distance(ind.begin(), w_it)-1;
            if(j<data.size()){
                vector<vector<double>> sub_vec(signal.begin()+i, signal.begin()+j+1);
                TimeSeries sub_signal(sub_vec);
                vector<int> id = get_id(sub_signal, w1, s1, delta, bsign);
                if (dict.count(id)){
                    dict[id].push_back(std::make_pair(i, j));
                    n++;
                }
                else{
                    vector<pair<int,int>> newlist {std::make_pair(i, j)};
                    dict.emplace(id, newlist);
                    n++;
                }
            }    
        }    
    }
    
    


}

/*
* Computes all relevant keys for a given query and collects the corresponding candidate solutions in the form of pairs of indices.
* @param query The query time series
* @return Returns candidate solutions in the form of pairs of indices.
*/
const std::vector<std::pair<int,int>> SubAnnDS::getAnswersIds1D(TimeSeries query){
    //this is only called during query
    vector<double> x = query.projectAxis(0);
    std::vector<int> key = get_query_frechet_id(x, w1,  s1, delta);
    std::vector<std::pair<int,int>> pair_ids, new_pair_ids; 
    for (int i = 0; i < 2; i++) {
        for (int j = key.size()-1; j < key.size(); j++) {
            std::vector<int> subkey(key.begin() + i, key.begin() + j + 1);
            new_pair_ids = dict[subkey];
            pair_ids.insert(pair_ids.end(), new_pair_ids.begin(), new_pair_ids.end());       
        }
    }
    //remove duplicates 
    sort(pair_ids.begin(), pair_ids.end(), lexcomparePairs);
    auto it = std::unique(pair_ids.begin(), pair_ids.end());
    pair_ids.erase(it, pair_ids.end());
    return pair_ids;
}


/*
* Range query for a given query.
* @param query The query time series
* @param R Distance threshold. Detected subcurves with distance at most R from the query will be returned.
* @return Returns detected subcurves in the form of pairs of indices. 
* A returned pair (i,j) indicates that there is a subcurve starting at 
* the (i-1)th edge and ending at the jth edge which is within distance R from the query
*/
std::vector<std::pair<int, int>> SubAnnDS::range_query(TimeSeries query, double R){
    std::vector<std::pair<int,int>> filtered_pair_ids;
    std::vector<std::pair<int,int>> pair_ids; 
    std::vector<TimeSeries> subqueries;    
    vector<double> x = query.projectAxis(0);
    std::vector<int> ind = signal.get_index();
    leq = false;
    //remove duplicates 
    pair_ids = getAnswersIds1D(query);
    std::map<std::pair<int, int>, bool> mask;
    
    #pragma omp parallel for 
    for (int j =0; j<pair_ids.size(); j++){
        std::pair<int,int> pair_id;

        auto pair_id0 = pair_ids[j];
        std::pair<int,int> pair_id1;
        TimeSeries starting_edge;
        TimeSeries ending_edge;
        // instead of pair_id.first-1, detect 
        std::vector<int> query_key = get_query_frechet_id(x, w1,  s1, delta);
        //TODO: EXTEND FOLLOWING TO MULTI-D
        
        // compute first key element of candidate        
        int ckey1 = get_point_id(signal[pair_id.first][0], w1, s1);
        // compute last key element of candidate
        int ckey2 = get_point_id(signal[pair_id.second][0], w1, s1);
        int i = pair_id0.first;
        int previous_key = ckey1;
        int new_key;
        //The following adds a "slack" of w1
        if(query_key[0]<=ckey1){
            while(i>=0){
                new_key = get_point_id(signal[i][0], w1, s1);
                if((new_key <= previous_key)&&(new_key >= query_key[0])){
                    i = i-1;
                    previous_key = new_key;
                }
                else{
                    break;
                }
            }
        }
        else{
            while(i>=0){
                new_key = get_point_id(signal[i][0], w1, s1);
                if((new_key >= previous_key)&&(new_key <= query_key[0])){
                    i = i-1;
                    previous_key = new_key;
                }
                else{
                    break;
                }
            }
        }
        
        pair_id1.first = i+1;
        i = pair_id0.second;
        previous_key = ckey2;
        if(query_key[query_key.size()-1]<=ckey2){
            while(i<signal.size()-1){
                new_key = get_point_id(signal[i][0], w1, s1);
                if((new_key <= previous_key)&&(new_key >= query_key[query_key.size()-1])){
                    i = i+1;
                    previous_key = new_key;
                }
                else{
                    break;
                }
            }
        }
        else{
             while(i<signal.size()-1){
                new_key = get_point_id(signal[i][0], w1, s1);
                if((new_key >= previous_key)&&(new_key <= query_key[query_key.size()-1])){
                    i = i+1;
                    previous_key = new_key;
                }
                else{
                    break;
                }
            }
        }
        pair_id1.second = i-1;

        for(int u = pair_id1.first; u<=pair_id0.first; u++){
            for(int z = pair_id0.second; z<=pair_id1.second; z++){
                if((window>0)&&(!leq)){
                    if (ind[z]-ind[u] != window){
                        break;
                    }
                }
                if((window>0)&&(leq)){
                    if (ind[z]-ind[u] > window){
                        break;
                    }
                }
                pair_id = make_pair(u,z);
                
                if(mask.count(pair_id)){
                    break;
                }
                    
                if(pair_id.first > 0){
                    starting_edge = signal.subseries(pair_id.first-1, pair_id.first);
                }
                else{
                    starting_edge = signal.subseries(pair_id.first, pair_id.first);
                        
                }
                if(pair_id.second<signal.size()-1){
                    ending_edge = signal.subseries(pair_id.second, pair_id.second+1);
                }
                else{
                    ending_edge = signal.subseries(pair_id.second, pair_id.second);
                }
                int dim = query.get_dim();
                double t;
                //Calculate x
                std::vector<double> x(dim);
                if(starting_edge.size()==2){
                    t = dot_product(vec_diff(starting_edge[0],starting_edge[1]), vec_diff(starting_edge[0],query[0]))/dot_product(vec_diff(starting_edge[0], starting_edge[1]), vec_diff(starting_edge[0], starting_edge[1]));
                    if (t >= -1 && t <= 0) {
                        for (int i = 0; i < dim; i++) {
                            x[i] = -t*starting_edge[1][i] + (1+t)*starting_edge[0][i];
                        }
                    }else if (norm(vec_diff(starting_edge[0], query[0])) < norm(vec_diff(starting_edge[1], query[0]))) {
                        x = starting_edge[0];
                    } else {
                        x = starting_edge[1];
                    }
                }
                else{
                    x = starting_edge[0];
                }
                
                //Calculate y
                std::vector<double> y(dim);
                if(ending_edge.size()==2){
                    t = dot_product(vec_diff(ending_edge[0],ending_edge[1]), vec_diff(ending_edge[0],query[0]))/dot_product(vec_diff(ending_edge[0], ending_edge[1]), vec_diff(ending_edge[0], ending_edge[1]));
                    if (t >= -1 && t <= 0) {
                        for (int i = 0; i < dim; i++) {
                            y[i] = -t*ending_edge[1][i] + (1+t)*ending_edge[0][i];
                        }
                    }else if (norm(vec_diff(ending_edge[0], query[0])) < norm(vec_diff(ending_edge[1], query[0]))) {
                        y = ending_edge[0];
                    } else {
                        y = ending_edge[1];
                    }
                }
                else{
                    y = ending_edge[0];
                }    
                //candidate curve
                auto start_it = signal.begin() + pair_id.first;
                auto end_it = signal.begin() + pair_id.second;
                // Create new vector with the desired subsequence
                vector<vector<double>> subseq(start_it, end_it+1);
                subseq.insert(subseq.begin(), vector<double>{x});
                subseq.insert(subseq.end(), vector<double>{y});
                
                if(Frechet::Continuous::less_than_or_equal(R, Curve(subseq), Curve(vector<vector<double>>(query.begin(), query.end())))){    
                
                    #pragma omp critical
                    {
                        filtered_pair_ids.push_back(pair_id);
                        mask[pair_id] = true;
                    }
                }
                
            }

        
        }    

    }
    //erase duplicates
    sort(filtered_pair_ids.begin(), filtered_pair_ids.end(), lexcomparePairs);
    auto it = std::unique(filtered_pair_ids.begin(), filtered_pair_ids.end());
    filtered_pair_ids.erase(it, filtered_pair_ids.end());
    return filtered_pair_ids;
    
}

/*
* Exact range query for a given query.
* @param query The query time series
* @param R Distance threshold. Detected subcurves with distance at most R from the query will be returned.
* @return Returns detected subcurves in the form of pairs of indices. 
* A returned pair (i,j) indicates that there is a subcurve starting at 
* the (i-1)th edge and ending at the jth edge which is within distance R from the query
*/
std::vector<std::pair<int, int>> SubAnnDS::exact_range_query(TimeSeries query, double R){
    std::pair<int,int> pair_id;
    std::vector<std::pair<int,int>> pair_ids; 
    std::vector<std::pair<int,int>> filtered_pair_ids;

    std::vector<TimeSeries> subqueries;    
    
    for (const auto& pair : dict) {
        for ( const auto x : pair.second){
           pair_ids.push_back(x);
       }
    }
    if (pair_ids.empty()){
        return std::vector<std::pair<int, int>>(1, make_pair(-1,-1));
    }
    std::vector<bool> mask(pair_ids.size(), false);
    #pragma omp parallel for 
    for (int j =0; j<pair_ids.size(); j++){
        std::pair<int,int> pair_id = pair_ids[j];
        TimeSeries starting_edge;
        TimeSeries ending_edge;
        // instead of pair_id.first-1, detect 
       
        //SETUP STARTING/ENDING EDGES

        if(pair_id.first > 0){
            starting_edge = signal.subseries(pair_id.first-1,pair_id.first);
        }
        else{
            starting_edge = signal.subseries(pair_id.first,pair_id.first);
            
        }
        if(pair_id.second<signal.size()-1){
            ending_edge = signal.subseries(pair_id.second,pair_id.second+1);
        }
        else{
            ending_edge = signal.subseries(pair_id.second,pair_id.second);
        }
        

        int dim = query.get_dim();
        double t;
        //Calculate x
        std::vector<double> x(dim);
        if(starting_edge.size()==2){
            t = dot_product(vec_diff(starting_edge[0],starting_edge[1]), vec_diff(starting_edge[0],query[0]))/dot_product(vec_diff(starting_edge[0], starting_edge[1]), vec_diff(starting_edge[0], starting_edge[1]));
            if (t >= -1 && t <= 0) {
                for (int i = 0; i < dim; i++) {
                    x[i] = -t*starting_edge[1][i] + (1+t)*starting_edge[0][i];
                }
            }else if (norm(vec_diff(starting_edge[0], query[0])) < norm(vec_diff(starting_edge[1], query[0]))) {
                x = starting_edge[0];
            } else {
                x = starting_edge[1];
            }
        }
        else{
            x = starting_edge[0];
        }
        
        //Calculate y
        std::vector<double> y(dim);
        if(ending_edge.size()==2){
            t = dot_product(vec_diff(ending_edge[0],ending_edge[1]), vec_diff(ending_edge[0],query[0]))/dot_product(vec_diff(ending_edge[0], ending_edge[1]), vec_diff(ending_edge[0], ending_edge[1]));
            if (t >= -1 && t <= 0) {
                for (int i = 0; i < dim; i++) {
                    y[i] = -t*ending_edge[1][i] + (1+t)*ending_edge[0][i];
                }
            }else if (norm(vec_diff(ending_edge[0], query[0])) < norm(vec_diff(ending_edge[1], query[0]))) {
                y = ending_edge[0];
            } else {
                y = ending_edge[1];
            }
        }
        else{
            y = ending_edge[0];
        }    
        //candidate curve
        auto start_it = signal.begin() + pair_id.first;
        auto end_it = signal.begin() + pair_id.second;
        // Create new vector with the desired subsequence
        vector<vector<double>> subseq(start_it, end_it+1);
        subseq.insert(subseq.begin(), vector<double>{x});
        subseq.insert(subseq.end(), vector<double>{y});
        if(Frechet::Continuous::less_than_or_equal(R, Curve(subseq), Curve(vector<vector<double>>(query.begin(), query.end())))){
            #pragma omp critical
            {
                filtered_pair_ids.push_back(pair_id);
            }
        }
        
    }
    return filtered_pair_ids;
    
}


/*
* Approximate nearest neighbor query for a given query.
* @param query The query time series
* @return Returns a subcurve in the form of pairs of indices. 
* A returned pair (i,j) indicates that there is a subcurve starting at 
* the (i-1)th edge and ending at the jth edge which is approximately the nearest to the query
*/
std::pair<int, int> SubAnnDS::ann_query(TimeSeries query){
    std::vector<std::pair<int,int>> filtered_pair_ids;
    std::vector<std::pair<int,int>> pair_ids; 
    std::vector<TimeSeries> subqueries;    
    vector<double> x = query.projectAxis(0);
    std::vector<int> ind = signal.get_index();
    leq = false;
    double min_dist = -1;
    std::pair<int, int> min_pair_id;
    //remove duplicates 
    pair_ids = getAnswersIds1D(query);
    std::map<std::pair<int, int>, bool> mask;
    #pragma omp parallel for 
    for (int j =0; j<pair_ids.size(); j++){
        std::pair<int,int> pair_id;

        auto pair_id0 = pair_ids[j];
        std::pair<int,int> pair_id1;
        TimeSeries starting_edge;
        TimeSeries ending_edge;
        std::vector<int> query_key = get_query_frechet_id(x, w1,  s1, delta);
        //TODO: EXTEND FOLLOWING TO MULTI-D
        
        // compute first key element of candidate        
        int ckey1 = get_point_id(signal[pair_id.first][0], w1, s1);
        // compute last key element of candidate
        int ckey2 = get_point_id(signal[pair_id.second][0], w1, s1);
        int i = pair_id0.first;
        int previous_key = ckey1;
        int new_key;
        //The following adds a "slack" of w1
        if(query_key[0]<=ckey1){
            while(i>=0){
                new_key = get_point_id(signal[i][0], w1, s1);
                if((new_key <= previous_key)&&(new_key >= query_key[0])){
                    i = i-1;
                    previous_key = new_key;
                }
                else{
                    break;
                    //TODO: break total
                }
            }
        }
        else{
            while(i>=0){
                new_key = get_point_id(signal[i][0], w1, s1);
                if((new_key >= previous_key)&&(new_key <= query_key[0])){
                    i = i-1;
                    previous_key = new_key;
                }
                else{
                    break;
                }
            }
        }
        
        pair_id1.first = i+1;
        i = pair_id0.second;
        previous_key = ckey2;
        if(query_key[query_key.size()-1]<=ckey2){
            while(i<signal.size()-1){
                new_key = get_point_id(signal[i][0], w1, s1);
                if((new_key <= previous_key)&&(new_key >= query_key[query_key.size()-1])){
                    i = i+1;
                    previous_key = new_key;
                }
                else{
                    break;
                }
            }
        }
        else{
             while(i<signal.size()-1){
                new_key = get_point_id(signal[i][0], w1, s1);
                if((new_key >= previous_key)&&(new_key <= query_key[query_key.size()-1])){
                    i = i+1;
                    previous_key = new_key;
                }
                else{
                    break;
                    //TODO: break total
                }
            }
        }
        pair_id1.second = i-1;

        for(int u = pair_id1.first; u<=pair_id0.first; u++){
            for(int z = pair_id0.second; z<=pair_id1.second; z++){
                pair_id = make_pair(u,z);
                if((window>0)&&(!leq)){
                    if(z<signal.size()-1){    
                        if(!((ind[z]-ind[u] <= window)&&(ind[z+1]-ind[u] > window))){///BUG
                                break;
                        }
                    }
                }
                if((window>0)&&(leq)){
                    if (ind[z]-ind[u] > window){//////// buuuugg
                        break;
                    }
                }
                
                if(mask.count(pair_id)){
                    break;
                }
                    
                if(pair_id.first > 0){
                    starting_edge = signal.subseries(pair_id.first-1, pair_id.first);
                }
                else{
                    starting_edge = signal.subseries(pair_id.first, pair_id.first);
                        
                }
                if(pair_id.second<signal.size()-1){
                    ending_edge = signal.subseries(pair_id.second, pair_id.second+1);
                }
                else{
                    ending_edge = signal.subseries(pair_id.second, pair_id.second);
                }
                int dim = query.get_dim();
                double t;
                //Calculate x
                std::vector<double> x(dim);
                if(starting_edge.size()==2){
                    t = dot_product(vec_diff(starting_edge[0],starting_edge[1]), vec_diff(starting_edge[0],query[0]))/dot_product(vec_diff(starting_edge[0], starting_edge[1]), vec_diff(starting_edge[0], starting_edge[1]));
                    if (t >= -1 && t <= 0) {
                        for (int i = 0; i < dim; i++) {
                            x[i] = -t*starting_edge[1][i] + (1+t)*starting_edge[0][i];
                        }
                    }else if (norm(vec_diff(starting_edge[0], query[0])) < norm(vec_diff(starting_edge[1], query[0]))) {
                        x = starting_edge[0];
                    } else {
                        x = starting_edge[1];
                    }
                }
                else{
                    x = starting_edge[0];
                }
                
                //Calculate y
                std::vector<double> y(dim);
                if(ending_edge.size()==2){
                    t = dot_product(vec_diff(ending_edge[0],ending_edge[1]), vec_diff(ending_edge[0],query[0]))/dot_product(vec_diff(ending_edge[0], ending_edge[1]), vec_diff(ending_edge[0], ending_edge[1]));
                    if (t >= -1 && t <= 0) {
                        for (int i = 0; i < dim; i++) {
                            y[i] = -t*ending_edge[1][i] + (1+t)*ending_edge[0][i];
                        }
                    }else if (norm(vec_diff(ending_edge[0], query[0])) < norm(vec_diff(ending_edge[1], query[0]))) {
                        y = ending_edge[0];
                    } else {
                        y = ending_edge[1];
                    }
                }
                else{
                    y = ending_edge[0];
                }    
                //candidate curve
                auto start_it = signal.begin() + pair_id.first;
                auto end_it = signal.begin() + pair_id.second;
                // Create new vector with the desired subsequence
                vector<vector<double>> subseq(start_it, end_it+1);
                subseq.insert(subseq.begin(), vector<double>{x});
                subseq.insert(subseq.end(), vector<double>{y});
       
                double dist = Frechet::Continuous::distance(Curve(subseq), Curve(vector<vector<double>>(query.begin(), query.end()))).value;
                #pragma omp critical
                {
                if ((min_dist == -1)||(dist < min_dist)){
                    min_dist = dist;
                    min_pair_id = pair_id;
                }
                }   
                
            }
        }    
    }
    //erase duplicates
    sort(filtered_pair_ids.begin(), filtered_pair_ids.end(), lexcomparePairs);
    auto it = std::unique(filtered_pair_ids.begin(), filtered_pair_ids.end());
    filtered_pair_ids.erase(it, filtered_pair_ids.end());
    return min_pair_id;
    
}



/*
* Exact nearest neighbor query for a given query.
* @param query The query time series
* @return Returns a subcurve in the form of pairs of indices. 
* A returned pair (i,j) indicates that there is a subcurve starting at 
* the (i-1)th edge and ending at the jth edge which is the nearest to the query
*/
std::pair<int, int> SubAnnDS::exact_nn_query(TimeSeries query){
    std::vector<std::pair<int,int>> filtered_pair_ids;
    std::vector<std::pair<int,int>> pair_ids; 
    std::vector<TimeSeries> subqueries;    
    vector<double> x = query.projectAxis(0);
    std::vector<int> ind = signal.get_index();
    leq = false;
    double min_dist = -1;
    std::pair<int, int> min_pair_id;
    //remove duplicates 
    for (const auto& pair : dict) {
        for ( const auto x : pair.second){
            pair_ids.push_back(x);
        }
    }
    std::map<std::pair<int, int>, bool> mask;
   
    #pragma omp parallel for 
    for (int j =0; j<pair_ids.size(); j++){
        std::pair<int,int> pair_id = pair_ids[j];
        
        TimeSeries starting_edge;
        TimeSeries ending_edge;
            
        if(pair_id.first > 0){
            starting_edge = signal.subseries(pair_id.first-1, pair_id.first);
        }
        else{
            starting_edge = signal.subseries(pair_id.first, pair_id.first);
                
        }
        if(pair_id.second<signal.size()-1){
            ending_edge = signal.subseries(pair_id.second, pair_id.second+1);
        }
        else{
            ending_edge = signal.subseries(pair_id.second, pair_id.second);
        }
        int dim = query.get_dim();
        double t;
        //Calculate x
        std::vector<double> x(dim);
        if(starting_edge.size()==2){
            t = dot_product(vec_diff(starting_edge[0],starting_edge[1]), vec_diff(starting_edge[0],query[0]))/dot_product(vec_diff(starting_edge[0], starting_edge[1]), vec_diff(starting_edge[0], starting_edge[1]));
            if (t >= -1 && t <= 0) {
                for (int i = 0; i < dim; i++) {
                    x[i] = -t*starting_edge[1][i] + (1+t)*starting_edge[0][i];
                }
            }else if (norm(vec_diff(starting_edge[0], query[0])) < norm(vec_diff(starting_edge[1], query[0]))) {
                x = starting_edge[0];
            } else {
                x = starting_edge[1];
            }
        }
        else{
            x = starting_edge[0];
        }
        
        //Calculate y
        std::vector<double> y(dim);
        if(ending_edge.size()==2){
            t = dot_product(vec_diff(ending_edge[0],ending_edge[1]), vec_diff(ending_edge[0],query[0]))/dot_product(vec_diff(ending_edge[0], ending_edge[1]), vec_diff(ending_edge[0], ending_edge[1]));
            if (t >= -1 && t <= 0) {
                for (int i = 0; i < dim; i++) {
                    y[i] = -t*ending_edge[1][i] + (1+t)*ending_edge[0][i];
                }
            }else if (norm(vec_diff(ending_edge[0], query[0])) < norm(vec_diff(ending_edge[1], query[0]))) {
                y = ending_edge[0];
            } else {
                y = ending_edge[1];
            }
        }
        else{
            y = ending_edge[0];
        }    
        //candidate curve
        auto start_it = signal.begin() + pair_id.first;
        auto end_it = signal.begin() + pair_id.second;
        // Create new vector with the desired subsequence
        vector<vector<double>> subseq(start_it, end_it+1);
        subseq.insert(subseq.begin(), vector<double>{x});
        subseq.insert(subseq.end(), vector<double>{y});
        
        double dist = Frechet::Continuous::distance(Curve(subseq), Curve(vector<vector<double>>(query.begin(), query.end()))).value;
        #pragma omp critical
        {
        if ((min_dist == -1)||(dist < min_dist)){
            min_dist = dist;
            min_pair_id = pair_id;
        }
        }   
        
    


    }
    //erase duplicates
    sort(filtered_pair_ids.begin(), filtered_pair_ids.end(), lexcomparePairs);
    auto it = std::unique(filtered_pair_ids.begin(), filtered_pair_ids.end());
    filtered_pair_ids.erase(it, filtered_pair_ids.end());
    return min_pair_id;
    
}

/*
*Outputs the contents of the dictionary 'dict' to the 'os' stream
*/
ostream& operator<<(ostream& os, const SubAnnDS& obj) {
    std::unordered_map<std::vector<int>, std::vector<std::pair<int,int>>, VectorHasher> dict = obj.get_dict();
    os<<"Printing elements of the dictionary"<<endl;
    for (const auto& pair : dict) {
        for (auto x: pair.first){
            os<<x;
        }    
        cout<< ": "; 
        for (auto x: pair.second){
            os<<"("<<x.first<<", "<<x.second<<")";
        } 
        os<< endl;
    }
    return os;
}    

/*
*Outputs basic statistics about the data structure
*/
void SubAnnDS::stats(){
    std::vector<int> n_buckets;
    for (auto p : dict) {
        n_buckets.push_back(p.second.size());
    }
    cout<<"Number of keys stored: "<<dict.size()<<endl;
    
    if(dict.size()>0){
        double sum = std::accumulate(n_buckets.begin(), n_buckets.end(), 0);
        double mean = sum / n_buckets.size();
        std::cout << "Total number of subsequences stored: " << sum << std::endl;

        std::cout << "Average number of subsequences per key: " << mean << std::endl;

        // // Compute maximum
        int max = *std::max_element(n_buckets.begin(), n_buckets.end());
        std::cout << "Maximum bucket size: " << max << std::endl;
        
    }    
}
