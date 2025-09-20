#ifndef CTSPC_H_
#define CTSPC_H_

#include <omp.h>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <cstring>
#include <algorithm>
#include <utility>
#include <cstdio>
#include <cmath>
#include <cstdint>
#include <climits>
#include <iostream>
#include <functional>
#include <unordered_map>

using namespace std;

#define MAXN 1000000000
#define MAXINT ((unsigned) 4294967295)// =2^32
#define MAXGID 100000000
#define MAXLINE 1024
#define MAXT 65000
#define MAXD 120

#define RANK_STATIC 0
#define RANK_LOCAL_STATIC 1 //LS
#define RANK_HOP_BETWEENNESS 2 //HB

typedef unsigned short tint;
typedef char dint;//指定标识符dint代表char类型

vector<int> score, _score;

class IntVal{
public:
    int x;
public:
    IntVal();
    IntVal(int _x);
    bool operator<(const IntVal &v) const;
};


class DV {
public:
    DV();
    DV(int id, double val);
    int id;
    double val;
    bool operator < (const DV& v) const;
};

class CD {
public:
    unsigned int cnt;
    int dis;
public:
    CD();
    CD(unsigned int cnt, int dis);
    bool operator<(const CD &v) const;
};


class TreeNode {
public:
    int id, f, h, rid, rsize, w;
    vector<int> nbr, cost, count; //of size <=w
    vector<int> anc; //of size h
    vector<dint> dis;
    vector<int> cnt;
    vector<int> ch;
};

class LabelEntry {
public:
    uint64_t v_d_c;
public:
    LabelEntry();
    LabelEntry(uint64_t _v);
    bool operator<(const LabelEntry &v) const;
};

class CoreTree {
public:
    static inline bool get_edge(char *line, int &a, int &b, int num_cnt = 2);
    static inline int get_num_cnt(string path);
    static void create_bin(string path,int rank_threads = 1,
            int rank_method = RANK_STATIC,  int rank_max_minutes = 1000000, int max_hops = 3, bool merge_equv = true);
    static void get_order(  vector<int> *con, int n, int *o, int method, int rank_threads, int rank_max_minutes, int max_hops);
    static void get_order_ls(vector<int> *con, int n, int *o, int rank_threads, int rank_max_minutes);
    static void get_order_hb(vector<int> *con, int n, int *o, int rank_threads, int rank_max_minutes, int n_hop);

public:
    string path;
    int **con, *dat, *deg, *nid, *oid, *f, *equ, *equ_c;
    bool *local;
    int n_org, n, n_core;
    long long m, m_core;
    unsigned MAXDIS, MAXMOV, MASK;
public:
    CoreTree(string path);
    ~CoreTree();
    void load_graph();

public:
    //vector<vector<int>> nbr, cost;
    vector<vector<int>> nbr, cost, count;
    vector<int> ord,rank;//ord:被max_w选中作为tree节点的点集 rank：没被选中的rank[i]=-1
    vector<vector<pair<int,CD>>> E;//存储 mde后的边信息
    vector<TreeNode> tree;
    double t;
    vector<pair<unsigned,unsigned>> *label;
    vector<pair<unsigned,unsigned>> *labelc;
    vector<pair<unsigned,unsigned>> *labeld;
    int max_w;

public:
    char can_update(int v, int dis, char *nowdis);
    void reduce(int max_w, int n_threads);
    void create_tree();

    void compute_tree_label();
    void compute_tree_label(int x, int rsize, vector<TreeNode*> &s);
    void compute_core_label(int n_threads);

    void decompose_tree(int max_w, int n_threads);
    void decompose_core(int max_w, int n_threads);

    void load_label_core(int max_w);
    void load_label_tree(int max_w);
    void load_label(int max_w);

    void save_label_core(int max_w);
    void save_label_tree(int max_w);

public:
    pair<int,CD> **core_con;
    pair<int,CD> *core_dat;
    int *core_deg;

    void save_tmp_graph(int max_w);
    void load_tmp_graph(int max_w);

public:
    tint *last_t;
    int *dis, * coun;
    tint nowt;

public:
    void init_query();
    int query(int u, int v);
    pair<int,int> query_by_core(int u, int v);
    pair<int,int> query_by_tree(int u, int v);
};

bool DV::operator < (const DV& v) const {
    if( val > v.val+1e-8 ) return true;
    if( val < v.val-1e-8 ) return false;
    return id > v.id;
}

DV::DV(int id, double val) {
    this->id = id; this->val = val;
}

DV::DV() {id = -1; val=-1;}


bool CD::operator < (const CD& v) const {
    if( cnt > v.cnt+1e-8 ) return true;
    if( cnt < v.cnt-1e-8 ) return false;
    return dis < v.dis;
}

CD::CD(unsigned int cnt, int dis) {
    this->cnt = cnt; this->dis = dis;
}

CD::CD() {cnt = 0; dis=-1;}



IntVal::IntVal() {x=-1;}
IntVal::IntVal(int x) {this->x=x;}
bool IntVal::operator<(const IntVal &v) const {if(score[x] == score[v.x]) return v.x<x; return score[x] < score[v.x];}

CoreTree::CoreTree(string path) {
    deg = NULL; dat = NULL; con = NULL;
    nid = NULL; oid = NULL; f = NULL;
    equ = NULL; equ_c=0;
    local = NULL;
    n_org = 0; n = 0; m = 0; n_core = 0; m_core = 0;
    MAXDIS = 0; MAXMOV = 0; MASK = 0; max_w = 0;
    label = NULL; 
    labelc = NULL; labeld = NULL;
    last_t = NULL; dis = NULL; nowt = 0; coun =NULL;
    core_con = NULL; core_dat = NULL;
    core_deg = NULL;
    this->path = path; t = omp_get_wtime();
}

CoreTree::~CoreTree() {
    if(deg) delete[] deg; if(dat) delete[] dat; if(con) delete[] con;
    if(nid) delete[] nid; if(oid) delete[] oid; if(f) delete[] f;
    if(equ) delete[] equ; if(equ_c) delete[] equ_c;
    if(local) delete[] local; 
    if(label) delete[] label; 
    if(labelc) delete[] labelc; if(labeld) delete[] labeld; 
    if(last_t) delete[] last_t; if(dis) delete[] dis;
    if(coun) delete[] coun;
    if(core_con) delete[] core_con; if(core_dat) delete[] core_dat;
    if(core_deg) delete[] core_deg; 
}

void CoreTree::init_query() {
    nowt = 0;
    last_t = new tint[n_core];
    memset( last_t, 0, sizeof(tint) * n_core );
    dis = new int[n_core];
    coun = new int[n_core];
}

pair<int,int> CoreTree::query_by_core(int u, int v) {
    ++nowt;
    
    if(nowt == MAXT) {
        memset(last_t, 0, sizeof(tint) * n_core);
        nowt=1;
    }
    
    if(rank[v]==-1 && rank[u]>=0) {int x=u; u=v; v=x;}
    if(rank[v]>=0 && rank[u]>=0)
        if(tree[v].rsize < tree[u].rsize) {int x=u; u=v; v=x;}
    if(rank[v]==-1 && rank[u]==-1)
        if(label[v].size() < label[u].size()) {int x=u; u=v; v=x;}
    int mind = MAXD;
    int mcount =0;
    if(rank[u] == -1) {
        for(auto l:label[u]) {
            last_t[l.first>>MAXMOV] = nowt; 
            dis[l.first>>MAXMOV] = l.first&MASK;
            coun[l.first>>MAXMOV] = l.second;
            }
    } else {
        TreeNode &tu = tree[u];
        TreeNode &r = tree[tu.rid];
        for(size_t i=0; i < tu.rsize; ++i) {
            int x = r.nbr[i], w = tu.dis[i], cou = tu.count[i];
            for(auto l:label[x]) {
                int y = l.first>>MAXMOV, nowd = (l.first&MASK)+w;
                unsigned nowc = l.second * cou* equ[y] ;
                if(nowd > mind) break;
                if(last_t[y] != nowt) {
                    last_t[y] = nowt; 
                    dis[y] = nowd;
                    coun[y] = nowc;
                }else if(nowd<dis[y]) {
                    dis[y] = nowd;
                    coun[y] = nowc;
                }else if(nowd==dis[y]) {
                    coun[y] += nowc;
                }
            }
        }
    }
    if(rank[v] == -1) {
        for(auto l:label[v]) {
            if(last_t[l.first>>MAXMOV] == nowt){
                int d = l.first&MASK + dis[l.first>>MAXMOV];
                if ( d < mind ){
                    mind = d; mcount = l.second * coun[l.first>>MAXMOV];
                }else if( d == mind){
                    mcount += l.second * coun[l.first>>MAXMOV];
                }
            }
        }
    } else {
        TreeNode &tv = tree[v];
        TreeNode &r = tree[tv.rid];
        for(size_t i = 0; i < tv.rsize; ++i ) {
            int x = r.nbr[i], w = tv.dis[i], cou = tv.count[i];
            for(auto l:label[x]) {
                int nowd = (l.first&MASK)+w, nowc = l.second * cou;
                if(nowd > mind) break;
                if(last_t[l.first>>MAXMOV] == nowt){
                    int d = nowd+dis[l.first>>MAXMOV];
                    if (d<mind){
                        mind = d;
                        mcount = nowc * coun[l.first>>MAXMOV];
                    } else if (d==mind){
                        mcount += nowc * coun[l.first>>MAXMOV];
                    }
                } 
            }
        }
    }
    return make_pair(mind,mcount);
}

pair<int,int> CoreTree::query_by_tree(int u, int v) {
    TreeNode &tu = tree[u], &tv = tree[v];
    int len = min(tu.h,tv.h);
    int d = INT_MAX;
    int num_c = 0;
    for(int i = tu.rsize, j = 0; i < len && tu.anc[j] == tv.anc[j]; ++i,++j){
        if (tu.dis[i]+tv.dis[i] < d){
            d = tu.dis[i]+tv.dis[i];
            num_c += tu.count[i]*tv.count[i];
        }else if(tu.dis[i]+tv.dis[i] == d){
            num_c += tu.count[i]*tv.count[i];
        }
    }
    return make_pair(d,num_c);
}


int CoreTree::query(int u, int v) {
    if( u == v ) return 1;
    int type = 0;
    u = nid[u]; v = nid[v];
    
    if(u < 0) {u = -u-1; type = 1;} else if(u >= MAXN) {u = u-MAXN; type = 2;}
    if(v < 0) {v = -v-1; type = 1;} else if(v >= MAXN) {v = v-MAXN; type = 2;}
    
    if(u == v) {
        if (type == 1){
            if (deg[u]==0){
                return 0;
            }else return equ_c[u];
        }else return 1;
    }
    
    //查询的vertex不合格
    if( u >= n || v >= n ) return 0;

    pair<int,int> result=query_by_core(u,v);
    int d = result.first;
    int c = result.second;
    //是同一棵树
    if(rank[u] >= 0 && rank[v] >= 0 && tree[u].rid == tree[v].rid){
        pair<int,int> tree_result=query_by_tree(u,v);
        if(tree_result.first <d){
            d = tree_result.first;
            c = tree_result.second;
        }else{
            c += tree_result.second;
        }
    }
    return c;
}

void CoreTree::save_tmp_graph(int max_w) {
    printf( "Saving Tmp Graph...\n" );
    char stw[16];
    sprintf(stw, "%d", max_w);
    FILE *fout = fopen( (path+"tmp-" + string(stw) + ".bin").c_str(), "wb" );

    fwrite(&n, sizeof(int), 1, fout);
    fwrite(&m_core, sizeof(long long), 1, fout);
    fwrite(rank.data(), sizeof(int), n, fout);
    vector<int> cored(n,0);
    for(int i = 0; i < n; ++i)
        if(rank[i] == -1) cored[i] = (int) E[i].size();
    fwrite(cored.data(), sizeof(int), n, fout);
    for(int i = 0; i < n; ++i)
        if(rank[i] == -1)
            fwrite(E[i].data(), sizeof(pair<int,CD>), E[i].size(), fout);
    fclose(fout);
    printf( "Tmp Graph Saved!\n" );
}

void CoreTree::load_tmp_graph(int max_w) {
    printf( "Loading Tmp Graph...\n" );
    char stw[16];
    sprintf(stw, "%d", max_w);

    FILE *fin = fopen( (path+"tmp-" + string(stw) + ".bin").c_str(), "rb" );
    fread(&n, sizeof(int), 1, fin);
    n_core = 0;
    fread(&m_core, sizeof(long long), 1, fin);
    rank.resize(n); 
    fread(rank.data(), sizeof(int), n, fin);
    for(int i = 0; i < n; ++i) if(rank[i] == -1) {++n_core; }
    printf( "n_core=%d\n", n_core );
    core_deg = new int[n]; core_dat = new pair<int,CD>[m_core]; core_con = new pair<int,CD>*[n];
    fread(core_deg, sizeof(int), n, fin);
    fread(core_dat, sizeof(pair<int,CD>), m_core, fin);
    long long p = 0;
    for(int i = 0; i < n; ++i)
        if(rank[i] ==-1) {
            core_con[i] = core_dat + p;
            p += core_deg[i];
        }
    printf( "m_core=%lld\n", p );
    fclose(fin);
    printf( "Tmp Graph Loaded!\n" );
}

void CoreTree::save_label_tree(int max_w) {
    printf( "Saving Tree Label...\n" );
    char stw[16];
    sprintf(stw, "%d", max_w);

    FILE *fout = fopen( (path+"label-tree-" + string(stw) + ".bin").c_str(), "wb" );
    fwrite(&n, sizeof(int), 1, fout);
    fwrite(rank.data(), sizeof(int), n, fout);
    for(int i = 0; i < n; ++i)
        if(rank[i] >= 0) {
            TreeNode &tn = tree[i];
            fwrite(&tn.rid, sizeof(int), 1, fout);
            fwrite(&tn.rsize, sizeof(int), 1, fout);
            fwrite(&tn.h, sizeof(int), 1, fout);
            fwrite(&tn.w, sizeof(int), 1, fout);
            fwrite(tn.nbr.data(), sizeof(int), tn.w, fout);
            fwrite(tn.anc.data(), sizeof(int), tn.h-tn.rsize, fout);
            fwrite(tn.dis.data(), sizeof(dint), tn.h, fout);
            fwrite(tn.cnt.data(), sizeof(int), tn.h, fout);
        }
    fclose(fout);
    printf( "Tree Label Saved!\n" );
}

void CoreTree::save_label_core(int max_w) {
    printf( "Saving Core Label...\n" );
    char stw[16];
    sprintf(stw, "%d", max_w);
    FILE *fout = fopen( (path+"label-core-" + string(stw) + ".bin").c_str(), "wb" );

    fwrite(&n, sizeof(int), 1, fout);
    //d
    for(int i = 0; i < n; ++i) {
        int len = (int) labeld[i].size();
        fwrite(&len, sizeof(int), 1, fout);
    }
    for(int i = 0; i < n; ++i)
        if(labeld[i].size() > 0)
            fwrite(labeld[i].data(), sizeof(pair<unsigned,unsigned>), labeld[i].size(), fout); 
    //c
    for(int i = 0; i < n; ++i) {
        int len = (int) labelc[i].size();
        fwrite(&len, sizeof(int), 1, fout);
    }
    for(int i = 0; i < n; ++i)
        if(labelc[i].size() > 0)
            fwrite(labelc[i].data(), sizeof(pair<unsigned,unsigned>), labelc[i].size(), fout); 

    fwrite(&MAXMOV, sizeof(int), 1, fout);
    fclose(fout);
    printf( "Core Label Saved!\n" );
}

void CoreTree::load_label_tree(int max_w) {
    printf( "Loading Tree Label...\n" );
    char stw[16];
    sprintf(stw, "%d", max_w);

    FILE *fin = fopen( (path+"label-tree-" + string(stw) + ".bin").c_str(), "rb" );
    fread(&n, sizeof(int), 1, fin);
    rank.resize(n);
    fread(rank.data(), sizeof(int), n, fin);
    n_core = 0;
    for(int i = 0; i < n; ++i) if(rank[i]==-1) ++n_core;
    tree.resize(n);
    for(int i = 0; i < n; ++i){
        if(rank[i] >= 0) {
            TreeNode &tn = tree[i];
            fread(&tn.rid, sizeof(int), 1, fin);
            fread(&tn.rsize, sizeof(int), 1, fin);
            fread(&tn.h, sizeof(int), 1, fin);
            fread(&tn.w, sizeof(int), 1, fin);
            tn.nbr.resize(tn.w); tn.anc.resize(tn.h); tn.dis.resize(tn.h);
            tn.count.resize(tn.h);
            fread(tn.nbr.data(), sizeof(int), tn.w, fin);
            fread(tn.anc.data(), sizeof(int), tn.h-tn.rsize, fin);
            fread(tn.dis.data(), sizeof(dint), tn.h, fin);
            fread(tn.count.data(), sizeof(int), tn.h, fin);
        }
    }
    fclose(fin);
}

void CoreTree::load_label_core(int max_w) {
    printf( "Loading Core Label...\n" );
    char stw[16];
    
    sprintf(stw, "%d", max_w);  

    FILE *fin = fopen( (path+"label-core-" + string(stw) + ".bin").c_str(), "rb" );
    fread(&n, sizeof(int), 1, fin);

    int *lend = new int[n];
    fread(lend, sizeof(int), n, fin);
    labeld = new vector<pair<unsigned,unsigned>>[n];
    
    label = new vector<pair<unsigned,unsigned>>[n];
    for( int i = 0; i < n; ++i ) {
        labeld[i].resize(lend[i]);
        fread(labeld[i].data(), sizeof(pair<unsigned,unsigned>), lend[i], fin);
        label[i].insert(label[i].end(),labeld[i].begin(),labeld[i].end());
    }
    int *lenc = new int[n];
    fread(lenc, sizeof(int), n, fin);
    labelc = new vector<pair<unsigned,unsigned>>[n];
    for( int i = 0; i < n; ++i ) {
        labelc[i].resize(lenc[i]);
        fread(labelc[i].data(), sizeof(pair<unsigned,unsigned>), lenc[i], fin);
        label[i].insert(label[i].end(),labelc[i].begin(),labelc[i].end());
        sort(label[i].begin(),label[i].end());
    }
    delete[] lend; delete[] lenc;
    
    fread(&MAXMOV, sizeof(unsigned), 1, fin);
    MAXDIS = 1<<MAXMOV; MASK = MAXDIS-1;
    fclose(fin);
    printf( "Core Label Loading Finished...\n" );
}

void CoreTree::load_label(int max_w) {
    this->max_w = max_w;
    load_label_tree(max_w);
    load_label_core(max_w);
    load_graph();
    init_query();
}

void CoreTree::decompose_tree(int max_w, int n_threads) {
    if(con == NULL) load_graph();
    reduce(max_w, n_threads);
    create_tree();
    compute_tree_label();
    save_label_tree(max_w);
}

void CoreTree::decompose_core(int max_w, int n_threads) {
    if(con == NULL) load_graph();
    load_tmp_graph(max_w);
    compute_core_label(n_threads);
    save_label_core(max_w);
}

char CoreTree::can_update(int v, int dis, char *nowdis) {
    int min=2;
    for(auto l:labeld[v]) {
        int d = l.first&MASK, h = l.first>>MAXMOV;
        if( nowdis[h] >= 0 ){
            if (nowdis[h] + d < dis) return -1;
            if (nowdis[h] + d == dis) min = 1;
        } 
    }
    return min;
}

void CoreTree::compute_core_label(int n_threads) {
    printf( "Computing Core Label...\n" );
    omp_set_num_threads(n_threads);
    if(n_core == 0) {printf("No core nodes!\n"); return;}
    MAXDIS = 2; MAXMOV = 1;
    while( MAXINT / (n_core * 2) >= MAXDIS ) {MAXDIS *= 2;++MAXMOV;}
    MASK = MAXDIS - 1;
    printf( "MAXDIS=%d\n", MAXDIS );
    //MAXDIS=536870912 = 2^29

    vector<int> *posd = new vector<int> [n];
    vector<int> *posc = new vector<int> [n];
    labelc = new vector<pair<unsigned,unsigned>>[n];
    labeld = new vector<pair<unsigned,unsigned>>[n];
    int p = 0;
    vector<int> vid;//进入core的 order
    vector<int> cid(n);//从原来的 --> core id
    for(int u = 0; u < n; ++u){
        if(rank[u]==-1 ) {
            //是core 
            vid.push_back(u);
            //进入core的 vertex original id
            cid[u] = p;
            labeld[u].push_back(make_pair( (((p++)<<MAXMOV) | 0) , 1));
            //移位操作自动转化成二进制
            posd[u].push_back(1);
            posc[u].push_back(0);
        } 
    }
    printf( "n_bc=%d,n_core=%d,n=%d\n", p, n_core, n );
    int dis = 1;
    for( long long cnt = 1; cnt && dis <= (int) MAXDIS; ++dis ) {
        cnt = 0;
        vector<pair<unsigned,unsigned>> *label_newc = new vector<pair<unsigned,unsigned>>[n];
        vector<pair<unsigned,unsigned>> *label_newd = new vector<pair<unsigned,unsigned>>[n];
        #pragma omp parallel
        {
            int pid = omp_get_thread_num(), np = omp_get_num_threads();
            long long local_cnt = 0;
            vector<char> used(n_core, 0);
            vector<int> candc;
            vector<int> candd;
            vector<int> candp;
            char *nowdis = new char[n_core];
            memset( nowdis, -1, sizeof(char) * n_core );
            unsigned *nowcount = new unsigned[n_core];
            memset( nowcount, 0, sizeof(unsigned) * n_core );
            #pragma omp for schedule(dynamic)
            for( int u = 0; u < n; ++u ) {
                if(rank[u] >= 0 ) continue;
                for(auto &l:labeld[u]) {
                    nowdis[l.first>>MAXMOV] = l.first&MASK;
                }
                candc.clear();
                candd.clear();
                candp.clear();
                for(int i = 0; i < core_deg[u]; ++i) {
                    int x = core_con[u][i].first, w = core_con[u][i].second.dis;
                    unsigned int wc = core_con[u][i].second.cnt;
                    if(w > dis) continue;
                    for( int j = (w==dis?0:posd[x][dis-w-1]); j < posd[x][dis-w]; ++j ) {
                        int v = labeld[x][j].first>>MAXMOV, vw=labeld[x][j].first &MASK, vc=labeld[x][j].second;       
                        if(vid[v] >= u) continue;
                        if (used[v] == 0){
                            used[v]=can_update(vid[v],dis,nowdis);
                            if (used[v] == 1){
                                candc.push_back(v);
                            }else if (used[v] == 2){
                                candd.push_back(v);
                            }else if (used[v] == -1){
                                candp.push_back(v);
                            }
                        }
                        if (used[v] >= 1){
                            nowcount[v] = vc*wc + nowcount[v];
                        }
                    }
                    for( int j = (w==dis?0:posc[x][dis-w-1]); j < posc[x][dis-w]; ++j ) {
                        int v = labelc[x][j].first>>MAXMOV, vw=labelc[x][j].first &MASK, vc=labelc[x][j].second;       
                        if(vid[v] >= u) continue;
                        if (used[v] == 0){
                            used[v]=can_update(vid[v],dis,nowdis);
                            if (used[v] == 1){
                                candc.push_back(v);
                            }else if (used[v] == 2){
                                candd.push_back(v);
                            }else if (used[v] == -1){
                                candp.push_back(v);
                            }
                        }
                        if (used[v] >= 1){
                            nowcount[v] = vc*wc + nowcount[v];
                        }
                    }
                }
                sort(candd.begin(), candd.end());
                sort(candc.begin(), candc.end());
                for(auto v:candd) used[v]=0;
                for(auto v:candc) used[v]=0;
                for(auto v:candp) used[v]=0;
                
                //size_t p = 0;
                for(auto v:candd) {
                    if(label_newd[u].size() > 100 && label_newd[u].size() == label_newd[u].capacity()) label_newd[u].reserve(label_newd[u].capacity() * 1.2);
                    label_newd[u].push_back(make_pair(  ( ( ( (unsigned)v) <<MAXMOV) | (unsigned) dis),nowcount[v] ) );
                    ++local_cnt;
                }
                for(auto v:candc) {
                    if(label_newc[u].size() > 100 && label_newc[u].size() == label_newc[u].capacity()) label_newc[u].reserve(label_newc[u].capacity() * 1.2);
                    label_newc[u].push_back(make_pair(  ( ( ( (unsigned)v) <<MAXMOV) | (unsigned) dis),nowcount[v] ) );
                    ++local_cnt;
                }
                
                for( int i = 0; i < (int) labeld[u].size(); ++i ) {
                    nowdis[(labeld[u][i].first>>MAXMOV)] = -1;
                }
                for( int i = 0; i < (int) label_newc[u].size(); ++i ) {
                    nowcount[(label_newc[u][i].first>>MAXMOV)] = 0;
                }
                for( int i = 0; i < (int) label_newd[u].size(); ++i ) {
                    nowcount[(label_newd[u][i].first>>MAXMOV)] = 0;
                }
            }

            if(pid==0) printf( "num_thread=%d,", np );
            #pragma omp critical
            {
                cnt += local_cnt;
            }
            delete[] nowdis;
            delete[] nowcount;
        }

        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic)
            for( int u = 0; u < n; ++u ) {
                if(rank[u] >= 0) continue;
                labeld[u].insert(labeld[u].end(), label_newd[u].begin(), label_newd[u].end());
                vector<pair<unsigned,unsigned>>().swap(label_newd[u]);
                labelc[u].insert(labelc[u].end(), label_newc[u].begin(), label_newc[u].end());
                vector<pair<unsigned,unsigned>>().swap(label_newc[u]);
                posd[u].push_back((int)labeld[u].size());
                posc[u].push_back((int)labelc[u].size());
            }
        }
        delete[] label_newc; delete[] label_newd;
        printf( "dis=%d,cnt=%lld,t=%0.3lf secs\n", dis, cnt,  omp_get_wtime()-t );
    }

    double tt = 0, max_label = 0;
    for( int i = 0; i < n_core; ++i ){
        if(rank[i] == -1) {
            if (local[i] == false){
                tt += labelc[i].size() ;
                tt += labeld[i].size() ;
            }
            max_label = max(max_label, (labeld[i].size()+ labelc[i].size())*1.0);
        }
    }
    printf( "Core label size=%0.3lfMB, Max core label size=%0.0lf, Avg core label size=%0.3lf, Time = %0.3lf sec\n",
            tt*8/(1024*1024.0), max_label, tt*0.125/n, omp_get_wtime()-t );
    
    delete[] posc; delete[] posd;
}

void CoreTree::compute_tree_label(int x, int rsize, vector<TreeNode*> &s) {
    s.push_back(&tree[x]);
    TreeNode &tn = tree[x];
    tn.dis.resize(tn.h);
    tn.cnt.resize(tn.h);
    int pos = 0;
    vector<int> p(tn.w);
    for(int j = 0; j < tn.w; ++j) {
        while(s[pos]->id != tn.nbr[j]) ++pos;
        p[j] = pos;
    }
    for(int i = 0; i < tn.h-1; ++i) {
        tn.dis[i] = -1;
        tn.cnt[i] = 0;
        int tv = (*s[i]).id;
        for(int j = 0; j < tn.w; ++j) {
            int w = tn.cost[j], k = p[j], nowdis = -1, nowcnt = 0,count = tn.count[j];
            int nv = tn.nbr[j];
            if(k<=i) {
                if(k==i){
                    nowdis=0;
                    nowcnt=1;
                }else if(i>=rsize) {
                    nowdis = s[i]->dis[k];
                    nowcnt = s[i]->dis[k];
                }
            }else if(k>=rsize) {
                nowdis = s[k]->dis[i],nowcnt  = s[k]->cnt[i];
            }
            if(nowdis>=0 && tv<=nv && (tn.dis[i]==-1 || nowdis+w<tn.dis[i])){
                tn.dis[i]=min(nowdis+w,MAXD);
                tn.cnt[i] = count * nowcnt;
            } else if(nowdis>=0 && tv<=nv && (tn.dis[i]==-1 || nowdis+w==tn.dis[i])){
                tn.cnt[i] += count * nowcnt;
            }
        }
    }
    tn.dis[tn.h-1] = 0;
    tn.cnt[tn.h-1] = 1;
    
    for(int &u:tree[x].ch)
        compute_tree_label(u, rsize, s);
    s.pop_back();
}

void CoreTree::compute_tree_label() {
    printf( "Computing Tree Label...\n" );
    vector<TreeNode*> s;
    for(int v=0; v<n; ++v)
        if(rank[v] >= 0 && tree[v].f == -1) {
            s.clear();
            for(int i=0; i<tree[v].rsize; ++i) s.push_back(&tree[tree[v].nbr[i]]);
            compute_tree_label(v, tree[v].rsize, s);
        }
    double t_size = 0;
    int maxdis = 0;
    for(int v=0; v<n; ++v) {
        if(rank[v] >= 0) {
            TreeNode &tn = tree[v];
            if (tn.ch.size()==0) continue;
            t_size += tree[v].dis.size() * 1.0 * (sizeof(int)+sizeof(int));
            vector<pair<int,CD>>().swap(E[v]);
            for(auto &d:tree[v].dis) maxdis = max(maxdis, (int)d);
        } else vector<pair<int,CD>>(E[v]).swap(E[v]);
    }
    printf( "Tree Label Computed, t=%0.3lf secs, maxdis=%d, tree label size=%0.3lf MB\n", omp_get_wtime() -t, maxdis, t_size/(1024.0*1024.0));
}

void CoreTree::create_tree() {
    printf( "Creating Tree...\n" );
    tree.resize(n);
    for(int u = 0; u < n; ++u) tree[u].id = u;
    vector<pair<int,int>> v_pair;
    int maxh = 0, cnt_root = 0, maxdep = 0, max_sub_tree = 1;
    vector<int> tcnt(n,0);
    
    double tw = 0;
    for(int i = (int) ord.size()-1; i >= 0; --i) {
        int x = ord[i];
        TreeNode &tn = tree[x];
        v_pair.clear();
        for(int j = 0; j < (int) nbr[x].size(); ++j) {
            int y = nbr[x][j];
            if(rank[y] == -1) v_pair.push_back(make_pair(n,j));
            else v_pair.push_back(make_pair(rank[y],j));
        }
        sort(v_pair.begin(),v_pair.end());
        reverse(v_pair.begin(), v_pair.end());
        
        int w = (int) nbr[x].size();
        tn.nbr.resize(w);
        tn.cost.resize(w);
        tn.count.resize(w);
        for(int j=0; j<w; ++j) {
            tn.nbr[j] = nbr[x][v_pair[j].second];
            tn.cost[j] = cost[x][v_pair[j].second];
            tn.count[j] = count[x][v_pair[j].second];
        }
        
        tn.w = w; 
        tn.id = x;
        tn.f = -1;
        for(auto &u:nbr[x])
            if(rank[u]!=-1 && (tn.f==-1 || rank[u] < rank[tn.f]))
                tn.f = u;
        if(tn.f == -1) {
            tn.h = tn.w + 1;
            ++cnt_root;//tree 的个数
            ++tcnt[x];//以x为根结点的tree node 的个数
            tn.rid = x;//root id
            tn.rsize = tn.w; // tree 对应的 bridge的size
            tn.anc.push_back(x);
        } else {
            tn.h = tree[tn.f].h+1;
            tree[tn.f].ch.push_back(x);
            tn.rid = tree[tn.f].rid;
            ++tcnt[tn.rid];
            max_sub_tree = max(max_sub_tree, tcnt[tn.rid]);
            tn.rsize = tree[tn.f].rsize;
            tn.anc = tree[tn.f].anc;
            tn.anc.push_back(x);
        }
        tw += tn.rsize;
        maxh = max(maxh, tn.h);
        maxdep = max(maxdep, (int)tn.anc.size());
    }
    printf( "Core tree constructed, maxh=%d, maxdep=%d, cnt_root=%d, max_stree=%d, avg_rsize=%0.3lf, t=%0.3lf secs\n",
          maxh, maxdep, cnt_root, max_sub_tree, tw/(n-n_core), omp_get_wtime()-t);
}

void CoreTree::reduce(int max_w, int n_threads) {
    omp_set_num_threads(n_threads);
    this->max_w = max_w;
    score.resize(n); _score.resize(n);
    vector<bool> changed(n,false);
    set<IntVal> q;
    nbr.resize(n); cost.resize(n); count.resize(n);
    rank.resize(n); fill(rank.begin(), rank.end(), -1);
    
    int r = 0;

    for(int i = 0; i < n; ++i) score[i] = _score[i] = deg[i];

    printf( "Initializing q..." );
    vector<bool> active(n,false);
    for(int u = 0; u < n; ++u){
        if(deg[u]<=max_w) {
            q.insert(IntVal(u));
            active[u] = true;
            }
        }
    printf( ", t=%0.3lf secs\nInitializing tmp...", omp_get_wtime()-t );
    E.resize(n);
    vector<unordered_map<int,CD>> tmp(n);
    for(int u = 0; u < n; ++u)
        for(int i = 0; i < deg[u]; ++i) {tmp[u][con[u][i]]={1,1};}
    printf( ", t=%0.3lf secs\nReducing Graph...\n", omp_get_wtime()-t );

    int cnt = 0;
    
    while(!q.empty()) {
        int x = q.begin()->x;
        while(changed[x]) {
            q.erase(x);
            score[x] = _score[x];
            q.insert(x);
            changed[x] = false;
            x = q.begin()->x;
        }
        if(score[x] > max_w) break;
        ord.push_back(x);
        q.erase(x);
        rank[x] = r++;
        //1. create Node v
        for (auto &e: tmp[x]) {
            if (rank[e.first] == -1) {
                E[x].emplace_back(e); 
            }
        }
        sort(E[x].begin(),E[x].end());
        for (auto & e: E[x]){
            nbr[x].push_back(e.first); count[x].push_back(e.second.cnt); cost[x].push_back(e.second.dis);
        }
        tmp[x].clear();
        //2. create clique
        for (const auto &s: E[x]) {
            int deg_inc = -1;
            for (const auto &item: E[x]) {
                if (s.first != item.first) {
                    auto it = tmp[s.first].find(item.first);
                    if (it == tmp[s.first].end()) {//new
                        tmp[s.first][item.first]={s.second.cnt * item.second.cnt, s.second.dis + item.second.dis};
                        deg_inc++;
                    }else {//old -->update
                        if (it->second.dis > s.second.dis + item.second.dis){
                            it->second = {s.second.cnt * item.second.cnt, s.second.dis + item.second.dis};
                        }else if (it->second.dis == s.second.dis + item.second.dis){
                            
                            it->second={tmp[s.first][item.first].cnt + s.second.cnt * item.second.cnt,tmp[s.first][item.first].dis};
                        }
                    }
                }
            }
            if (deg_inc != 0) {
                _score[s.first] = _score[s.first] + deg_inc;
                if(active[s.first] ==false && (rank[s.first] == -1)){
                    score[s.first] = _score[s.first];
                }
                if(_score[s.first] >= max_w * 2) {active[s.first] = false; q.erase(s.first); continue;}
                changed[s.first] = true;
            }
        }
        if((++cnt) * score[x] > 1000000) {
            printf( "%d nodes reduced, score[x]=%d, remaining size=%0.3lf%% t=%0.3lf secs\n",
                    r, (n-r)*100.0/n, score[x], omp_get_wtime()-t);
            cnt = 0;
        }
    }
    
    printf( "Reordering edges...\n" );
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for (int u = 0; u < n; ++u) {
            if (rank[u] != -1) continue;
            vector<pair<int, CD>> tmp_edges;
            for (auto &e: tmp[u]) {
                if (rank[e.first] != -1) continue;
                tmp_edges.emplace_back(e);
            }
            sort(tmp_edges.begin(), tmp_edges.end());
            swap(tmp_edges, E[u]);
        }
    }
    n_core = 0;
    m_core = 0;
    for(int u = 0; u < n; ++u)
        if( rank[u] == -1 ) {
            ++n_core;
            m_core += (int) E[u].size();
        }
    printf( "Reducing finished, t=%0.3lf secs, n_core=%d,m_core=%lld,node_rate=%0.3lf,edge_rate=%0.3lf\n",
            omp_get_wtime()-t, n_core, m_core, n_core*1.0/n, m_core*1.0/m );
}

void CoreTree::load_graph() {
    printf( "loading graph: %s\n", path.c_str() );
    long long p = 0;
    FILE* fin = fopen( (path+"graph-dis.bin").c_str(), "rb" );
    fread( &n, sizeof(int), 1, fin );
    fread( &m, sizeof(long long), 1, fin );
    deg = new int[n]; dat = new int[m]; con = new int*[n]; nid = new int[n]; oid = new int[n];
    local = new bool[n];
    equ = new int[n], equ_c = new int[n];
    fread( deg, sizeof(int), n, fin );
    fread( dat, sizeof(int), m, fin );
    fread( nid, sizeof(int), n, fin );
    fread( oid, sizeof(int), n, fin );
    fread( equ, sizeof(int), n, fin );
    fread( equ_c, sizeof(int), n, fin );
    fread( local, sizeof(bool), n, fin );
    fclose(fin);
    for( int i = 0; i < n; ++i ) {con[i] = dat+p; p+= deg[i];}
    int nown = n-1; while(nown>=0 && deg[nown] == 0) --nown;
    nown += 1; if( nown < 2 ) nown = 2;
    n_org = n; n = nown;
    printf( "graph loaded, n_org=%d, n=%d, m=%lld\n", n_org, n, m);
}

bool CoreTree::get_edge(char *line, int &a, int &b, int num_cnt) {
    if( !isdigit(line[0]) ) return false;
    vector<char*> v_num;
    int len = (int) strlen(line);
    for( int i = 0; i < len; ++i )
        if( !isdigit(line[i]) && line[i] != '.') line[i] = '\0';
        else if(i == 0 || !line[i-1]) v_num.push_back(line+i);
    if( (int) v_num.size() != num_cnt ) return false;
    sscanf( v_num[0], "%d", &a );
    sscanf( v_num[1], "%d", &b );
    return true;
}

int CoreTree::get_num_cnt(string path) {
    FILE *fin = fopen( (path + "graph.txt").c_str(), "r" );
    char line[MAXLINE];
    int cnt = 0, min_cnt = 100;

    while( fgets( line, MAXLINE, fin ) && cnt < 10 ) {
        if( !isdigit(line[0]) ) continue;
        vector<char*> v_num;
        int len = (int) strlen(line);
        for( int i = 0; i < len; ++i )
            if( !isdigit(line[i]) && line[i] != '.' ) line[i] = '\0';
            else if(i == 0 || !line[i-1]) v_num.push_back(line+i);
        if( (int) v_num.size() < 2 ) continue;
        min_cnt = min(min_cnt, (int) v_num.size());
        ++cnt;
    }
    fclose( fin );
    return min_cnt;
}

void CoreTree::get_order(vector<int> *con, int n, int *o, int method, int rank_threads, int rank_max_minutes, int max_hops) {
    printf( "method=%d\n", method );
    if( method == RANK_STATIC ) {
        printf( "Ranking Method = RANK_STATIC\n" );
        DV *f = new DV[n];
        for( int i = 0; i < n; ++i )
            f[i].id = i, f[i].val = con[i].size() * 1.0;
        sort(f, f + n);
        for(int i = 0; i < n; ++i) {
            o[i] = f[i].id;
        }
        delete[] f;
    } else if( method == RANK_LOCAL_STATIC ) {
        get_order_ls(con, n, o, rank_threads, rank_max_minutes);
    } else if( method == RANK_HOP_BETWEENNESS ) {
        get_order_hb(con, n, o, rank_threads, rank_max_minutes, max_hops);
    }
}

void CoreTree::get_order_ls(vector<int> *con, int n, int *o,  int rank_threads, int rank_max_minutes) {
    printf( "Ranking Method = RANK_LOCAL_STATIC\n" );

    int *ord = new int[n];
    for( int i = 0; i < n; ++i ) ord[i] = i;
    for(int i = n-1; i > 0; --i ) {
        int p = rand()%i;
        int tmp = ord[i]; ord[i] = ord[p]; ord[p] = tmp;
    }

    DV *f = new DV[n];

    for( int i = 0; i < n; ++i ) {
        f[i].id = i;
        f[i].val = con[i].size() == 0 ? 0 : 0.1;
    }

    int *cnt = new int[n];
    memset(cnt, 0, sizeof(int) * n);
    int s_used = 0;

    omp_set_num_threads(rank_threads);
    double t =  omp_get_wtime();
    #pragma omp parallel
    {
        int pid = omp_get_thread_num(), np = omp_get_num_threads();
        if( pid == 0 ) printf( "rank_threads = %d, rank_max_minutes = %d\n", np, rank_max_minutes );

        int *loc_r = new int[n];
        memset(loc_r, 0, sizeof(int) * n);

        vector<int> loc_pair;

        bool finished = false;

        #pragma omp for schedule(dynamic)
        for( int p = 0; p < n; ++p ) {
            if( finished ) continue;
            if( p % 100000 == 0 ) printf( "[%d] t=%0.3lf sec\n", p, omp_get_wtime()-t );
            int u = ord[p];
            loc_pair.clear();
            for( int i = 0; i < (int) con[u].size(); ++i ) {
                int v = con[u][i];
                for( int j = 0; j < (int) con[v].size(); ++j ) {
                    int w = con[v][j];
                    if( w >= u ) break;
                    if( loc_r[w] == 0 ) loc_pair.push_back(w);
                    ++loc_r[w];
                }
            }
            for( int i = 0; i < (int) con[u].size(); ++i ) {
                int v = con[u][i];
                double s_val = 0.0;
                for( int j = 0; j < (int) con[v].size(); ++j ) {
                    int w = con[v][j];
                    if( w >= u ) break;
                    double val = 1.0/loc_r[w];
                    s_val += val;
                }

                #pragma omp critical
                ++cnt[v], f[v].val += s_val;
            }

            #pragma omp critical
            ++s_used;

            for( int i = 0; i < (int) loc_pair.size(); ++i ) loc_r[loc_pair[i]] = 0;

            if( omp_get_wtime()-t >  rank_max_minutes * 60.0 ) finished = true;
        }

        delete[] loc_r;
    }

    printf( "Ranking Finished, ranking time=%0.3lf sec, sampling rate=%0.4lf%%\n", omp_get_wtime()-t, s_used*100.0/n );
    for( int i = 0; i < n; ++i )
        if(cnt[i]>0) f[i].val = f[i].val * (con[i].size()*1.0/cnt[i]);
    sort(f, f + n);
    for(int i = 0; i < n; ++i) {
        o[i] = f[i].id;
        if( i<10 ) printf( "id=%d,val=%0.3lf,cnt=%d,size=%d\n", f[i].id, f[i].val, cnt[o[i]], (int)con[o[i]].size());
    }

    delete[] f; delete[] ord; delete[] cnt;
}

void CoreTree::get_order_hb(vector<int> *con, int n, int *o, int rank_threads, int rank_max_minutes, int n_hop) {
    printf( "Ranking Method = RANK_HOP_BETWEENNESS\n" );

    int *ord = new int[n];
    for( int i = 0; i < n; ++i ) ord[i] = i;
    for(int i = n-1; i > 0; --i ) {
        int p = rand()%i;
        int tmp = ord[i]; ord[i] = ord[p]; ord[p] = tmp;
    }

    DV *f = new DV[n];

    for( int i = 0; i < n; ++i ) {
        f[i].id = i;
        f[i].val = con[i].size() == 0 ? 0 : 0.1;
    }

    int *cnt = new int[n];
    memset(cnt, 0, sizeof(int) * n);
    bool *used = new bool[n];
    memset(used, 0, sizeof(bool) * n);
    int s_used = 0;

    omp_set_num_threads(rank_threads);
    double t =  omp_get_wtime();
    #pragma omp parallel
    {
        int pid = omp_get_thread_num(), np = omp_get_num_threads();
        if( pid == 0 ) printf( "rank_threads = %d, rank_max_minutes = %d, max_hop=%d\n", np, rank_max_minutes, n_hop );

        double *loc_n = new double[n];
        double *loc_b = new double[n];
        char *loc_d = new char[n];
        memset(loc_d, -1, sizeof(char) * n);

        vector<int> loc_v;

        bool finished = false;

        #pragma omp for schedule(dynamic)
        for( int p = 0; p < n; ++p ) {
            if( finished ) continue;
            if( p % 100000 == 0 ) printf( "Stage 1 Sampling: [%d] t=%0.3lf sec\n", p, omp_get_wtime()-t );
            int s = ord[p];
            loc_v.clear(); loc_v.push_back(s);
            loc_d[s] = 0; loc_n[s] = 1;
            int pre_p = 0, now_p = 1;
            for( int h = 1; h <= n_hop; ++h ) {
                for( int i = pre_p; i < now_p; ++i ) {
                    int u = loc_v[i];
                    for( int j = 0; j < (int) con[u].size(); ++j ) {
                        int v = con[u][j];
                        if(loc_d[v] == -1) {
                            loc_n[v] = 0; loc_b[v] = 0;
                            loc_d[v] = h; loc_v.push_back(v);
                        }
                        if(loc_d[v] == loc_d[u] + 1) loc_n[v] += loc_n[u];
                    }
                }
                pre_p = now_p; now_p = (int) loc_v.size();
            }

            for( int i = (int) loc_v.size() - 1; i > 0; --i )
                if( loc_d[loc_v[i]] < n_hop ) {
                    int u = loc_v[i];
                    for( int j = 0; j < (int) con[u].size(); ++j ) {
                        int v = con[u][j];
                        if( loc_d[v] != loc_d[u] + 1 ) continue;
                        loc_b[u] += loc_n[u]/loc_n[v] * (1+loc_b[v]);
                    }

                    #pragma omp critical
                    f[u].val += loc_b[u]*1.0*loc_d[u], ++cnt[u];
            }

            for( int i = 0; i < (int) loc_v.size(); ++i ) loc_d[loc_v[i]] = -1;

            #pragma omp critical
            ++s_used, used[s] = true;

            if( omp_get_wtime()-t >  rank_max_minutes * 60.0 * 0.8 ) finished = true;
        }

        delete[] loc_n; delete[] loc_b; delete[] loc_d;
    }

    printf( "Stage 1 Sampling Finished!\n" );

    int *total = new int[n];
    memcpy( total, cnt, sizeof(int) * n );
    int s_cnt = s_used;

    #pragma omp parallel
    {
        vector<int> loc_v;
        bool finished = false;

        bool *loc_u = new bool[n];
        memset(loc_u, 0, sizeof(bool) * n);

        #pragma omp for schedule(dynamic)
        for( int p = 0; p < n; ++p ) {
            if( finished || used[ord[p]] ) continue;
            if( p % 100000 == 0 ) printf( "Stage 2 Counting: [%d] t=%0.3lf sec\n", p, omp_get_wtime()-t );
            int s = ord[p];
            loc_v.clear(); loc_v.push_back(s);
            loc_u[s] = true;
            int pre_p = 0, now_p = 1;
            for( int h = 1; h < n_hop; ++h ) {
                for( int i = pre_p; i < now_p; ++i ) {
                    int u = loc_v[i];
                    for( int j = 0; j < (int) con[u].size(); ++j ) {
                        int v = con[u][j];
                        if( !loc_u[v] ) {
                            loc_u[v] = true; loc_v.push_back(v);
                        }
                    }
                }
                pre_p = now_p; now_p = (int) loc_v.size();
            }

            #pragma omp critical
            for( int i = 1; i < (int) loc_v.size(); ++i )
                ++total[loc_v[i]];

            #pragma omp critical
            ++s_cnt;

            if( omp_get_wtime()-t >  rank_max_minutes * 60.0) finished = true;
            for( int i = 0; i < (int)loc_v.size(); ++i ) loc_u[loc_v[i]] = false;
        }
        delete[] loc_u;
    }

    printf( "Ranking Finished, ranking time=%0.3lf sec, sampling rate=%0.4lf%%, counting rate=%0.4lf%%\n",
            omp_get_wtime()-t, s_used*100.0/n, s_cnt*100.0/n );
    for( int i = 0; i < n; ++i )
        if(cnt[i]>0) f[i].val = f[i].val * (total[i]*1.0/cnt[i]);
    sort(f, f + n);
    for(int i = 0; i < n; ++i) {
        o[i] = f[i].id;
        if( i<10 ) printf( "id=%d,val=%0.3lf,total=%d,cnt=%d,deg=%d\n", f[i].id, f[i].val, total[o[i]], cnt[o[i]], (int)con[o[i]].size());
    }

    delete[] f; delete[] ord; delete[] cnt; delete[] total; delete[] used;
}


void CoreTree::create_bin(string path, int rank_threads ,int rank_method, int rank_max_minutes, int max_hops, bool merge_equv ) {
    FILE *fin = fopen( (path + "graph.txt").c_str(), "r" );
    char line[MAXLINE];
    int n = 0, a, b, num_cnt = get_num_cnt(path);
    vector< pair<int,int> > el;
    long long cnt = 0, m = 0;
    while( fgets( line, MAXLINE, fin ) ) {
        if( !get_edge(line, a, b, num_cnt) ) continue;
        if( a < 0 || b < 0 || a == b ) continue;
        el.push_back(make_pair(a, b));
        n = max(max(n, a+1), b+1);
        if( (++cnt) % (long long) 10000000 == 0 ) printf( "%lld lines finished\n", cnt );
    }
    fclose( fin );

    double nowt = omp_get_wtime();
    vector<int> *con = new vector<int>[n];
    printf( "Deduplicating...\n" );

    for(size_t i = 0; i < el.size(); ++i) {
        con[el[i].first].push_back(el[i].second);
        con[el[i].second].push_back(el[i].first);
    }
    for( int i = 0; i < n; ++i )
        if( con[i].size() > 0 ){
            sort( con[i].begin(), con[i].end() );
            int p = 1;
            for( int j = 1; j < (int) con[i].size(); ++j )
                if( con[i][j-1] != con[i][j] ) con[i][p++] = con[i][j];
            con[i].resize( p ); m += p;
        }

    long long *f1 = new long long[n];
    memset( f1, 0, sizeof(long long) * n );

    long long *f2 = new long long[n];
    memset( f2, 0, sizeof(long long) * n );

    
    if( !merge_equv ) {
        for( int i = 0; i < n; ++i ) f1[i] = i, f2[i] = i;
    } else {
        //printf( "Merging...\n" );
        long long s = 0;
        long long *nows = new long long[m+n+1];
        int *nowt = new int[m+n+1];
        memset( nowt, 0, sizeof(int) * (m+n+1) );

        for( int v = 0; v < n; ++v ){
            for( int i = 0; i < (int) con[v].size(); ++i ) {
                int u = con[v][i];
                if( nowt[f1[u]] != (v+1) ) {
                    ++s;
                    nows[f1[u]] = s;
                    nowt[f1[u]] = (v+1);
                    f1[u] = s;
                } else f1[u] = nows[f1[u]];
            }
        }
       
        for( int v = 0; v < n; ++v )
            if( nowt[f1[v]] != -1 ) {
                nows[f1[v]] = v;
                nowt[f1[v]] = -1;
                f1[v] = v;
            } else f1[v] = nows[f1[v]];

        s = 0;
        memset( nowt, 0, sizeof(int) * (m+n+1) );
        for( int v = 0; v < n; ++v )
            for( int i = 0; i <= (int) con[v].size(); ++i ) {
                int u = (i == (int) con[v].size()) ? v : con[v][i];
                if( nowt[f2[u]] != (v+1) ) {
                    ++s;
                    nows[f2[u]] = s;
                    nowt[f2[u]] = (v+1);
                    f2[u] = s;
                } else f2[u] = nows[f2[u]];
            }

        for( int v = 0; v < n; ++v )
            if( nowt[f2[v]] != -1 ) {
                nows[f2[v]] = v;
                nowt[f2[v]] = -1;
                f2[v] = v;
            } else f2[v] = nows[f2[v]];
                
        delete[] nows; delete[] nowt;

        long long cnt1_n = 0, cnt1_m = 0, cnt2_n = 0, cnt2_m = 0;
        for( int i = 0; i < n; ++i ) {
            if( f1[i] != i ) {
                ++cnt1_n;
                cnt1_m += (int) con[i].size();
            }
            if( f2[i] != i ) {
                ++cnt2_n;
                cnt2_m += (int) con[i].size();
            }
        }

        m = 0;
        for( int i = 0; i < n; ++i ) {
            if( f1[i] != i || f2[i] != i ) {con[i].clear(); continue;}
            int p = 0;
            for( int j = 0; j < (int) con[i].size(); ++j ) {
                int v = con[i][j];
                if( f1[v] == v && f2[v] == v ) con[i][p++] = v;
            }
            con[i].resize(p); m += p;
        }
        printf( "cnt1_n = %lld, cnt1_m = %lld, cnt2_n = %lld, cnt2_m = %lld, m = %lld\n", cnt1_n, cnt1_m, cnt2_n, cnt2_m, m );
    }



    printf( "Reordering...\n" );
    int *f = new int[n];
    get_order(con, n, f, rank_method, rank_threads, rank_max_minutes, max_hops);
    int *oid = new int[n], *nid = new int[n], *equ = new int[n], *equ_c = new int[n];

    for( int i = 0; i < n; ++i )
        oid[i] = f[i], nid[f[i]] = i, equ[i]=1, equ_c[i]=0;
   
    for( int i = 0; i < n; ++i ) {
        if(f1[i] != i) {
            nid[i] = -nid[f1[i]]-1;
            equ[nid[f1[i]]]++;
            equ_c[nid[f1[i]]] = con[f1[i]].size();
            }
        if(f2[i] != i){
            nid[i] = MAXN + nid[f2[i]];
            equ[nid[f2[i]]]++;
            equ_c[nid[f2[i]]] =1;
        } 
    }
    
    printf( "Creating adjacency list...\n" );
    int *dat = new int[m], *deg = new int[n], **adj = new int *[n];
    
    long long pos = 0;
    for( int i = 0; i < n; ++i ) {
        adj[i] = dat + pos;
        pos += (int) con[oid[i]].size();

    }
    memset( deg, 0, sizeof(int) * n );
    
    //deg=更新后的degree 对应 order之后的vertex 
    //++++++++dat对应的 相邻的点 位置
    for( int i = 0; i < n; ++i ) {
        int ii = oid[i];
        for( int p = 0; p < (int) con[ii].size(); ++p ) {
            int jj = con[ii][p];
            int j = nid[jj];
            adj[j][deg[j]++] = i;
        }
    }  
    bool *local = new bool[n];
    memset( local, false, sizeof(bool) * n );

    for (int u = 0; u < n; ++u) {
        if (deg[u] == 0) continue;
        int cnt = 0;
        for (int j =0;j<deg[u];j++){
            if (u > adj[u][j]) ++cnt;
        }
        if (deg[u] == cnt) local[u] = true;
    }
    
    nowt = omp_get_wtime() - nowt;
    printf( "Creating Bin Time = %0.3lf secs\n", nowt );

    printf( "Saving binary...\n" );
    FILE *fout = fopen( (path + "graph-dis.bin").c_str(), "wb" );


    fwrite( &n, sizeof(int), 1, fout );
    fwrite( &m, sizeof(long long), 1, fout );
    fwrite( deg, sizeof(int), n, fout );
    fwrite( dat, sizeof(int), m, fout );
    fwrite( nid, sizeof(int), n, fout );
    fwrite( oid, sizeof(int), n, fout );
    fwrite( equ, sizeof(int), n, fout );
    fwrite( equ_c, sizeof(int), n, fout );
    fwrite( local, sizeof(bool), n, fout );
    fclose( fout );
    
    printf( "Created txt file, n = %d, m = %lld\n", n, m );
    
    delete[] adj; delete[] deg; delete[] dat; delete[] f; delete[] con; delete[] oid; delete[] nid;
    delete[] f1; delete[] f2; delete[] equ; delete[] equ_c; delete[] local;
}

#endif /* CTSPC_H_ */
s
