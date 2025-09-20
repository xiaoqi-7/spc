#include "ctspc.h"
#include <fstream>

void make_query(string path, int n_pairs, string filename, int max_w){
    CoreTree ct(path);
    ct.load_graph();
    ct.load_tmp_graph(max_w);
    int n = ct.n_org;
    vector<int> rank = ct.rank;
    int *nid = ct.nid;
    int u=0,v=0;
    ofstream outfile_cores;
    outfile_cores.open(("qd/"+filename+"_cores"+to_string(max_w) +".txt").c_str());
    
    int type = 0;
    outfile_cores<<n_pairs<<"\n";
    int k =0;
    for(int i=0 ; k< n_pairs ; i++){
        int ou = rand() % n, ov = rand() % n;
        u = nid[ou];
        v = nid[ov];
        if(u < 0) {u = -u-1;} else if(u >= MAXN) {u = u-MAXN;}
        if(v < 0) {v = -v-1;} else if(v >= MAXN) {v = v-MAXN;}
        if(rank[u] == -1 && rank[v] == -1 && ou != ov){
            outfile_cores<<ou<<" "<<ov<<"\n";
            k=k+1;
        }
    }
    outfile_cores.close();
    cout<<"core finished"<<endl;
    if(max_w==0) return;
    ofstream outfile_ct;
    outfile_ct.open(("qd/"+filename+"_ct"+to_string(max_w)+".txt").c_str());
    
    k =0;
    for(int i=0 ; k< n_pairs ; i++){
        int ou = rand() % n, ov = rand() % n;
        u = nid[ou];
        v = nid[ov];
        if(u < 0) {u = -u-1;} else if(u >= MAXN) {u = u-MAXN;}
        if(v < 0) {v = -v-1;} else if(v >= MAXN) {v = v-MAXN;}
        if(rank[u] == -1 && rank[v] != -1 && ou != ov){
            outfile_ct<<ou<<" "<<ov<<"\n";
            k=k+1;
        }
    }
    outfile_ct.close();
    cout<<"ct finished"<<endl;
    
    ofstream outfile_dtt;
    outfile_dtt.open(("qd/"+filename+"_dtt"+to_string(max_w)+".txt").c_str());
    
    k =0;
    for(int i=0 ; k< n_pairs ; i++){
        int ou = rand() % n, ov = rand() % n;
        u = nid[ou];
        v = nid[ov];
        if(u < 0) {u = -u-1;} else if(u >= MAXN) {u = u-MAXN;}
        if(v < 0) {v = -v-1;} else if(v >= MAXN) {v = v-MAXN;}
        if(rank[u] != -1 && rank[v] != -1 && ov != ou){
            outfile_dtt<<ou<<" "<<ov<<"\n";
            k=k+1;
        }
    }
    outfile_dtt.close();
    cout<<"dtt finished"<<endl;

    ofstream outfile_stt;
    outfile_stt.open(("qd/"+filename+"_stt"+to_string(max_w)+".txt").c_str());
    
    k =0;
    for(int i=0 ; k< n_pairs ; i++){
        int ou = rand() % n, ov = rand() % n;
        u = nid[ou];
        v = nid[ov];
        if(u < 0) {u = -u-1;} else if(u >= MAXN) {u = u-MAXN;}
        if(v < 0) {v = -v-1;} else if(v >= MAXN) {v = v-MAXN;}
        if(rank[u] != -1 && rank[v] != -1 && ov != ou){
            outfile_stt<<ou<<" "<<ov<<"\n";
            k=k+1;
        }
    }
    outfile_stt.close();
    cout<<"stt finished"<<endl;
}


void query_count(string path, int n_pairs, string filename, int max_w){
    CoreTree ct(path);
    ct.load_label(max_w);
    FILE* file_cores = fopen(("qd/"+filename+"_cores"+to_string(max_w)+".txt").c_str(), "r");
    uint32_t num_queries = 0;
    fscanf(file_cores, "%d" , &num_queries);
    vector<std::pair<uint32_t, uint32_t>> queries_cores;
    for (uint32_t q = 0; q < num_queries; ++q) {
        uint32_t v1, v2;
        fscanf(file_cores, "%d"  "%d" , &v1, &v2);
        queries_cores.push_back({v1, v2});
    }
    fclose(file_cores);
    double t = omp_get_wtime();
    for(auto query:queries_cores){
        uint32_t v1 = query.first;
        uint32_t v2 = query.second;
        int count = ct.query(v1,v2);
    }
    t = omp_get_wtime() - t;
    
    printf("time_cores = %0.3lfns ", t*1000000.0/num_queries);
    cout<<endl;
    if(max_w==0) return;
    FILE* file_ct = fopen(("qd/"+filename+"_ct"+to_string(max_w)+".txt").c_str(), "r");
    num_queries = 0;
    fscanf(file_ct, "%d" , &num_queries);
    vector<std::pair<uint32_t, uint32_t>> queries_ct;
    for (uint32_t q = 0; q < num_queries; ++q) {
        uint32_t v1, v2;
        fscanf(file_ct, "%d"  "%d" , &v1, &v2);
        queries_ct.push_back({v1, v2});
    }
    fclose(file_ct);
    t = omp_get_wtime();
    for(auto query:queries_ct){
        uint32_t v1 = query.first;
        uint32_t v2 = query.second;
        int dis = ct.query(v1,v2);
    }
    t = omp_get_wtime() - t;
    printf("time_ct = %0.3lfns ", t*1000000.0/num_queries);
    cout<<endl;

    FILE* file_dtt = fopen(("qd/"+filename+"_dtt"+to_string(max_w)+".txt").c_str(), "r");
    num_queries = 0;
    fscanf(file_dtt, "%d" , &num_queries);
    vector<std::pair<uint32_t, uint32_t>> queries_dtt;
    for (uint32_t q = 0; q < num_queries; ++q) {
        uint32_t v1, v2;
        fscanf(file_dtt, "%d"  "%d" , &v1, &v2);
        queries_dtt.push_back({v1, v2});
    }
    fclose(file_dtt);
    t = omp_get_wtime();
    for(auto query:queries_dtt){
        uint32_t v1 = query.first;
        uint32_t v2 = query.second;
        int dis = ct.query(v1,v2);
    }
    t = omp_get_wtime() - t;
    printf("time_dtt = %0.3lfns ", t*1000000.0/num_queries);
    cout<<endl;

    FILE* file_stt = fopen(("qd/"+filename+"_stt"+to_string(max_w)+".txt").c_str(), "r");
    num_queries = 0;
    fscanf(file_stt, "%d" , &num_queries);
    vector<std::pair<uint32_t, uint32_t>> queries_stt;
    for (uint32_t q = 0; q < num_queries; ++q) {
        uint32_t v1, v2;
        fscanf(file_stt, "%d"  "%d" , &v1, &v2);
        queries_stt.push_back({v1, v2});
    }
    fclose(file_stt);
    t = omp_get_wtime();
    for(auto query:queries_stt){
        uint32_t v1 = query.first;
        uint32_t v2 = query.second;
        int dis = ct.query(v1,v2);
    }
    t = omp_get_wtime() - t;
    printf("time_stt = %0.3lfns ", t*1000000.0/num_queries);
    cout<<endl;
}

int main(int argc, char *argv[]) {
    setvbuf(stdout, NULL, _IONBF, 0);
    setvbuf(stderr, NULL, _IONBF, 0);
    if( argc > 1 ) {
        if(strcmp( argv[1], "txt-to-bin" ) == 0){
            CoreTree::create_bin(
                argv[2],                                /*path*/
                argc>3?atoi(argv[3]):1,                 /*rank_threads*/
                argc>4?atoi(argv[4]):RANK_STATIC,       /*rank_method*/
                argc>5?atoi(argv[5]):60,                /*rank_max_minutes*/
                argc>6?atoi(argv[6]):3                  /*max_hops*/
            );
        }else if(strcmp(argv[1], "decompose_tree") == 0) {
            CoreTree ct(argv[2]);                       /*path*/
            ct.decompose_tree(atoi(argv[3]),             /*max_w*/
                argc>4?atoi(argv[4]):1                  /*n_threads*/
            );  
            
        } else if(strcmp(argv[1], "decompose_core") == 0) {
            CoreTree ct(argv[2]);                       /*path*/
            ct.decompose_core(atoi(argv[3]),            /*max_w*/
                argc>4?atoi(argv[4]):1                  /*n_threads*/
            );
        } else if(strcmp(argv[1], "decompose_bt") == 0) {
            CoreTree ct(argv[2]);                       /*path*/
            ct.decompose_tree(atoi(argv[3]),            /*max_w*/
                argc>4?atoi(argv[4]):1                  /*n_threads*/
            );
            ct.save_tmp_graph(atoi(argv[3]));
        }else if (strcmp(argv[1], "make_queries") == 0){
            make_query(argv[2],                         /*path*/
                atoi(argv[3]),                          /*n_pairs*/
                argv[4],                                /*filename*/
                atoi(argv[5]));                         /*max_w*/
        }else if (strcmp(argv[1], "query_spc") == 0){
            query_count(argv[2],                        /*path*/
                atoi(argv[3]),                          /*n_pairs*/
                argv[4],                                /*filename*/
                atoi(argv[5]));                         /*max_w*/
        }
    }

    return 0;
}
