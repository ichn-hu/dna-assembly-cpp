#include "approx_match.cpp"
#include "genome.cpp"

struct Patch {
    int l, s;
    int dis, end;
    Patch(int l = 0, int s = 0, int dis = 0, int end = 0)
        : l(l)
        , s(s)
        , dis(dis)
        , end(end)
    {
    }
    void save(FILE* f)
    {
        fprintf(f, "%d %d %d %d\n", l, s, dis, end);
    }
    bool read(FILE* f)
    {
        int d = fscanf(f, "%d%d%d%d", &l, &s, &dis, &end);
        return d == 4;
    }
};

struct RepairFactory {
    vector<Genome> longReads;
    vector<Genome> shortReads;
    vector<Patch> patches;
    Config cfg;
    int lid, sid;
    RepairFactory(Config cfg)
        : cfg(cfg)
    {
        longReads = readFasta(cfg.longPath);
        auto s1 = readFasta(cfg.shortPath1);
        auto s2 = readFasta(cfg.shortPath2);
        for (unsigned int i = 0; i < s1.size(); ++i) {
            shortReads.push_back(s1[i]);
            shortReads.push_back(s2[i]);
        }
        puts("");
        for (sid = 0; sid < (int)shortReads.size(); sid += 2) {
            findPatch(&shortReads[sid]);
            sid += 1;
            findPatch(&shortReads[sid]);
            printf("\b\b\b\b\b\b\b\b\b\b\b");
            printf("%d/%d", sid, (int)shortReads.size());
            if (sid % 100 == 0) {
                exportPatches();
                printf("\n export up to %d\n", sid);
            }
            fflush(stdout);
            // break;
        }
        exportPatches();
    }
    void findPatch(Genome* s)
    {
        lid = 0;
        #pragma opm parallel for
        for (int i = 0; i < (int)longReads.size(); ++i) {
            findPatch(&longReads[i], s, lid);
            ++lid;
            // if (lid % 100 == 0) {
            //     printf("%d/%d\n", lid, (int)longReads.size());
            // }
        }
    }
    void findPatch(Genome* l, Genome* s, int lid)
    {
        auto S = l->data;
        auto T = s->data;
        int n = S.length(), m = T.length();
        int e[1000 + 1][100 + 1];
        memset(e, 0x3f, sizeof e);
        for (int i = 0; i <= n; ++i)
            e[0][i] = 0;
        for (int i = 0; i <= m; ++i)
            e[i][0] = i;
        for (int i = 1; i <= m; ++i) {
            for (int j = 1; j <= n; ++j) {
                if (e[i - 1][j] + 1 < e[i][j]) {
                    e[i][j] = e[i - 1][j] + 1;
                }
                if (e[i][j - 1] + 1 < e[i][j]) {
                    e[i][j] = e[i][j - 1] + 1;
                }
                if (T[i - 1] == S[j - 1]) {
                    if (e[i - 1][j - 1] < e[i][j]) {
                        e[i][j] = e[i - 1][j - 1];
                    }
                } else {
                    if (e[i - 1][j - 1] + 1 < e[i][j]) {
                        e[i][j] = e[i - 1][j - 1] + 1;
                    }
                }
            }
        }
        for (int i = 0; i <= n; ++i) {
            if (e[m][i] < (int)T.length() * 0.22) {
                patches.push_back(Patch(sid, lid, e[m][i], i));
            }
        }
    }

    void exportPatches()
    {
        auto f = fopen(cfg.patchesPath.c_str(), "w");
        for (auto&& p : patches) {
            p.save(f);
        }
    }
    void importPatches()
    {
        patches.clear();
        Patch t;
        auto f = fopen(cfg.patchesPath.c_str(), "r");
        while (t.read(f)) {
            patches.push_back(t);
        }
    }
};

main()
{
    auto cfg = Config("data4");
    auto fac = RepairFactory(cfg);
}