#include "config.h"
#include <algorithm>
#include <cstring>
#include <map>
#include <string>
#include <vector>
using namespace std;

char comp(char c)
{
    switch (c) {
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    }
}

int c2i(char c)
{
    switch (c) {
    case 'A':
        return 0;
    case 'T':
        return 1;
    case 'G':
        return 2;
    case 'C':
        return 3;
    }
}

char i2c(int i)
{
    switch (i) {
    case 0:
        return 'A';
    case 1:
        return 'T';
    case 2:
        return 'G';
    case 3:
        return 'C';
    }
}

// O(n^2) to compute the minimal edit distance between two strings
int _minEditDist(string a, string b)
{
    int n = a.length();
    int m = b.length();
    int f[n + 1][m + 1];
    memset(f, 0x3f, sizeof f);
    f[0][0] = 0;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; ++j) {
            f[i][j] = min(f[i][j], f[i - 1][j] + 1);
            f[i][j] = min(f[i][j], f[i][j - 1] + 1);
            if (a[i - 1] == b[j - 1]) {
                f[i][j] = min(f[i][j], f[i - 1][j - 1]);
            } else {
                f[i][j] = min(f[i][j], f[i - 1][j - 1] + 1);
            }
        }
    }
    return f[n][m];
}

class Genome {
public:
    string data;
    Genome() {}
    Genome(string t) { this->data.assign(t); }
    Genome(string::iterator b, string::iterator e) { this->data.assign(b, e); }
    Genome reverse()
    {
        Genome t(this->data);
        std::reverse(t.data.begin(), t.data.end());
        return t;
    }
    Genome compliment()
    {
        int l = this->data.length();
        char* s = new char[l + 1];
        memcpy(s, this->data.c_str(), l);
        for (int i = 0; i < l; ++i) {
            s[i] = comp(s[i]);
        }
        s[l] = 0;
        Genome t(s);
        return t;
    }
    int minEditDist(Genome t)
    {
        return _minEditDist(this->data, t.data);
    }
    bool operator<(const Genome rhs) const
    {
        return data < rhs.data;
    }
    bool operator==(const Genome rhs) const
    {
        return data == rhs.data;
    }
    vector<Genome> kmers(int k)
    {
        vector<Genome> ret;
        if (data.length() < k) {
            return ret;
        }
        for (unsigned int i = 0; i < data.length() - k; ++i) {
            ret.push_back(Genome(data.substr(i, k)));
        }
        return ret;
    }
    char lastChar()
    {
        return data[data.length() - 1];
    }
    char firstChar()
    {
        return data[0];
    }
    Genome firstKm1mer()
    {
        return Genome(data.substr(0, data.length() - 1));
    }
    Genome lastKm1mer()
    {
        return Genome(data.substr(1, data.length()));
    }
    vector<char> getRawData()
    {
        vector<char> ret;
        for (auto&& c : data) {
            ret.push_back(c);
        }
        return ret;
    }
};

vector<Genome> readFasta(string path)
{
    printf("reading from %s\n", path.c_str());
    auto f = fopen(path.c_str(), "r");
    char s[233333];
    vector<Genome> ret;
    while (fscanf(f, "%*s\n%s\n", s) != EOF) {
        auto g = Genome(s);
        ret.push_back(g);
        // printf("%s\n%s\n%s\n%s\n\n", g.data.c_str(), g.reverse().data.c_str(),
        // g.reverse().compliment().data.c_str(),
        // g.compliment().reverse().data.c_str());
    }
    printf("%d genomes read from %s\n", (int)ret.size(), path.c_str());
    return ret;
}

struct DBGNode;

struct DBGEdge {
    DBGNode* node;
    int cvg;
    DBGEdge(DBGNode* node = nullptr, int cvg = 0)
        : node(node)
        , cvg(cvg)
    {
    }
};

struct DBGNode : Genome {
    DBGEdge *from[4], *to[4];
    int id;
    bool visited;
    void init()
    {
        for (int i = 0; i < 4; ++i) {
            from[i] = nullptr;
            to[i] = nullptr;
        }
        visited = false;
    }
    DBGNode(Genome g)
    {
        data.assign(g.data);
        init();
    }
    int numTo()
    {
        int n = 0;
        for (int i = 0; i < 4; ++i)
            if (to[i] != nullptr)
                ++n;
        return n;
    }
    int numFrom()
    {
        int n = 0;
        for (int i = 0; i < 4; ++i) {
            if (from[i] != nullptr) {
                ++n;
            }
        }
        return n;
    }
};

struct DeBrujinGraph {
    vector<DBGNode*> nodes;
    map<Genome, int> toId;

    vector<int> extLen;

    DBGNode* getNode(Genome mer)
    {
        auto has = toId.find(mer);
        if (has == toId.end()) {
            auto newNode = new DBGNode(mer);
            newNode->id = nodes.size();
            toId[mer] = newNode->id;
            nodes.push_back(newNode);
            return newNode;
        }
        return nodes[has->second];
    }
    void addEdge(DBGEdge*& edges, DBGNode* node)
    {
        if (edges == nullptr) {
            edges = new DBGEdge(node, 0);
        }
        ++edges->cvg;
    }
    DeBrujinGraph(vector<Genome>& genomes)
    {
        printf("Building DBG... ");
        fflush(stdout);
        for (auto&& genome : genomes) {
            DBGNode *prev = nullptr, *curr = nullptr;
            for (auto&& kmer : genome.kmers(cfg.k)) {
                prev = getNode(kmer.firstKm1mer());
                curr = getNode(kmer.lastKm1mer());
                addEdge(prev->to[c2i(curr->lastChar())], curr);
                addEdge(curr->from[c2i(prev->firstChar())], prev);
            }
        }
        printf("done\n");
        fflush(stdout);
    }
    void analysis()
    {
        printf("Analysing...\n");
        int numHead = 0;
        int totNode;
        // int cntSucc[5] = { 0 };
        // int cntPrev[5] = { 0 };
        map<pair<int, int>, int> cnt;
        for (int i = 0; i <= 4; ++i) {
            for (int j = 0; j <= 4; ++j)
                cnt[make_pair(i, j)] = 0;
        }
        for (auto&& node : nodes) {
            // ++cntPrev[node->numFrom()];
            // ++cntSucc[node->numTo()];
            cnt[make_pair(node->numFrom(), node->numTo())]++;
            // printf("%s %s\n", node.second->data.c_str(), node.first.data.c_str());
        }
        printf("From/To: \n");
        for (int i = 0; i <= 4; ++i) {
            for (int j = 0; j <= 4; ++j) {
                printf("%7d", cnt[make_pair(i, j)]);
            }
            printf("\n");
        }
        // printf("%d/%d to --0: %6d --1: %6d --2: %6d --3: %6d --4: %6d\n", numHead, totNode, cntSucc[0], cntSucc[1], cntSucc[2], cntSucc[3], cntSucc[4]);
    }
    void removeBubble()
    {
    }
    int dfsExtLen(int currId)
    {
        if (extLen[currId] != 0) {
            return extLen[currId];
        }
        auto curr = nodes[currId];
        extLen[currId] = 1;
        int maxSucc = 0;
        for (int i = 0; i < 4; ++i) {
            auto to = curr->to[i];
            if (to != nullptr) {
                maxSucc = max(maxSucc, dfsExtLen(to->node->id));
            }
        }
        extLen[currId] += maxSucc;
        return extLen[currId];
    }
    void calcExtLen()
    {
        extLen.resize(nodes.size());
        for (auto&& i : extLen)
            i = 0;
        int maxLen = 0;
        for (auto i = 0u; i < nodes.size(); ++i) {
            maxLen = max(maxLen, dfsExtLen(i));
        }
        printf("Extend length: %d\n", maxLen);
        fflush(stdout);
    }
    void followMaxExtLen(DBGNode* curr, vector<char>& path)
    {
        if (extLen[curr->id] != 1) {
            for (int i = 0; i < 4; ++i) {
                if (curr->to[i] != nullptr) {
                    auto next = curr->to[i]->node;
                    if (extLen[next->id] == extLen[curr->id] - 1) {
                        path.push_back(i2c(i));
                        followMaxExtLen(next, path);
                        break;
                    }
                }
            }
        }
    }
    void output()
    {
        auto f = fopen(cfg.resultPath.c_str(), "w");
        printf("Output to %s\n", cfg.resultPath.c_str());
        fflush(stdout);
        int tot = 0, src = 0;
        calcExtLen();
        for (auto&& node : nodes) {
            ++tot;
            if (node->numFrom() == 0) {
                ++src;
                vector<char> path = node->getRawData();
                followMaxExtLen(node, path);
                printf("Outputing %d/%d: len = %d, extLen = %d\n", src, tot, (int)path.size(), extLen[node->id]);
                fflush(stdout);
                fprintf(f, ">%d/%d\n", src, tot);
                for (auto&& c : path) {
                    fprintf(f, "%c", c);
                }
                fprintf(f, "\n");
            }
        }
    }
};

void extendGenomes(vector<Genome>& genomes)
{
    auto tmp = vector<Genome>();
    for (auto&& g : genomes) {
        tmp.push_back(g.reverse());
    }
    genomes.reserve(genomes.size() + tmp.size());
    genomes.insert(genomes.end(), tmp.begin(), tmp.end());
    tmp.clear();
    for (auto&& g : genomes) {
        tmp.push_back(g.compliment());
    }
    genomes.reserve(genomes.size() + tmp.size());
    genomes.insert(genomes.end(), tmp.begin(), tmp.end());
}

main()
{
    auto genomes = readFasta(cfg.shortPath1);
    extendGenomes(genomes);
    auto graph = new DeBrujinGraph(genomes);
    graph->analysis();
    // graph->output();
    return 0;
}
