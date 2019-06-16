#include "config.h"
#include <algorithm>
#include <cstring>
#include <map>
#include <set>
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

void writeFasta(vector<Genome>& genomes, string path, int threshold)
{
    auto f = fopen(cfg.resultPath.c_str(), "w");
    printf("Outputing %d genomes to %s...", (int)genomes.size(), cfg.resultPath.c_str());
    fflush(stdout);
    int minLen = 0x3f3f3f3f, maxLen = 0, totLen = 0;
    for (int i = 0; i < (int)genomes.size(); ++i) {
        int len = (int)genomes[i].data.length();
        if (len < threshold) {
            continue;
        }
        minLen = min(minLen, len);
        maxLen = max(maxLen, len);
        totLen += len;
        fprintf(f, ">%d: len = %d\n", i, len);
        fprintf(f, "%s\n", genomes[i].data.c_str());
    }
    printf("done, average length %d, max %d, min %d\n", totLen / (int)genomes.size(), maxLen, minLen);
    fflush(stdout);
}

