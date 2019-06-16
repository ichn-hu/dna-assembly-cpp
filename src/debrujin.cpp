#include "genome.cpp"

struct DBGNode;

struct DBGEdge {
    DBGNode* node;
    int cvg;
    int leadToLoop;
    bool visited;
    bool removed;
    DBGEdge(DBGNode* node = nullptr, int cvg = 0)
        : node(node)
        , cvg(cvg)
    {
        leadToLoop = 0;
        visited = false;
        removed = false;
    }
};

struct DBGNode : Genome {
    DBGEdge *from[4], *to[4];
    int id;
    bool visited;
    bool removed;
    void init()
    {
        for (int i = 0; i < 4; ++i) {
            from[i] = nullptr;
            to[i] = nullptr;
        }
        visited = false;
        removed = false;
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
    Config cfg;

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
    DeBrujinGraph(vector<Genome>& genomes, Config cfg)
        : cfg(cfg)
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
    void analysisNodeInOut()
    {
        printf("Analysing node in out degree...\n");
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
            if (node->removed)
                continue;
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

    vector<multiset<int>> branches;
    vector<int> numBranches;
    vector<bool> visited;
    void dfsAnalysisBranch(DBGNode* u)
    {
        if (visited[u->id]) {
            return;
        }
        visited[u->id] = true;
        if (u->numTo() == 0) {
            branches[u->id].insert(u->id);
            numBranches[u->id] = 1;
            return;
        }
        for (int i = 0; i < 4; ++i) {
            auto v = u->to[i];
            if (v != nullptr) {
                dfsAnalysisBranch(v->node);
                numBranches[u->id] += numBranches[v->node->id];
                for (auto&& t : branches[v->node->id]) {
                    branches[u->id].insert(t);
                }
            }
        }
    }
    void analysisBranch()
    {
        printf("Analyising branches ...\n");
        branches.resize(nodes.size());
        fill(branches.begin(), branches.end(), multiset<int>());

        numBranches.resize(nodes.size());
        fill(numBranches.begin(), numBranches.end(), 0);

        visited.resize(nodes.size());

        for (auto&& u : nodes) {
            if (u->numFrom() == 0) {
                fill(visited.begin(), visited.end(), false);
                dfsAnalysisBranch(u);
                printf("Num branches for %d: %d\n", u->id, numBranches[u->id]); //(int)branches[u->id].size());
                for (auto&& t : branches[u->id]) {
                    printf("%d ", t);
                }
                printf("\n");
                fflush(stdout);
            }
        }
    }

    vector<int> timeStamp;
    // vector<bool> onStack;
    int visClk;
    int dfsFindLoop(DBGNode* u)
    {
        timeStamp[u->id] = ++visClk;
        // onStack[u->id] = true;
        int uLow = timeStamp[u->id];
        int numTo = u->numTo();
        for (int i = 0; i < 4; ++i) {
            auto v = u->to[i];
            if (v != nullptr) {
                if (timeStamp[v->node->id] == 0) {
                    int vLow = dfsFindLoop(v->node);
                    if (vLow < uLow) {
                        uLow = vLow;
                        if (numTo > 1) {
                            v->leadToLoop = vLow;
                        }
                    }
                } else {
                    uLow = min(uLow, timeStamp[v->node->id]);
                }
            }
        }
        return uLow;
    }
    void findLoop()
    {
        timeStamp.resize(nodes.size());
        fill(timeStamp.begin(), timeStamp.end(), 0);
        visClk = 0;
        for (auto&& u : nodes) {
            if (u->numFrom() == 0) {
                printf("Finding loop for node %d ... ", u->id);
                visClk = 0;
                fill(timeStamp.begin(), timeStamp.end(), 0);
                dfsFindLoop(u);
                printf(" done\n");
                fflush(stdout);
            }
        }
    }
    void analysisLoop()
    {
        printf("Analysing loop... ");
        fflush(stdout);
        findLoop();
        int loopCnt = 0;
        for (auto&& u : nodes) {
            for (int i = 0; i < 4; ++i) {
                auto v = u->to[i];
                if (v != nullptr && v->leadToLoop != 0) {
                    ++loopCnt;
                }
            }
        }
        printf("%d loops found\n", loopCnt);
        /*
        ➜  dna-assembly-cpp git:(master) ✗ make 1
        g++ src/genome.cpp -DDATA1 -g -o build/ass.exe
        build/ass.exe
        reading from ./data/data1/short_1.fasta
        8500 genomes read from ./data/data1/short_1.fasta
        Building DBG... done
        Analysing node in out degree...
        From/To:
            0     76      0      0      0
            76 168832    192      0      0
            0    192      0      0      0
            0      0      0      0      0
            0      0      0      0      0
        Analysing loop... 85 loops found
        */
        // Weird... Why 85 loops? It should be 192 loops...
    }
    void analysis()
    {
        analysisNodeInOut();
        // analysisLoop();
    }
    // if visited path has smaller length, then it must be a bubble, otherwise it might be a loop
    void walkThroughBubble(DBGNode* u, vector<pair<DBGNode*, int>>& path)
    {
        while (u->numTo() == 1) {
            int i;
            for (i = 0; i < 4; ++i)
                if (u->to[i] != nullptr) {
                    break;
                }
            path.push_back(make_pair(u, i));
            u = u->to[i]->node;
            if (u->numFrom() != 1)
                break;
            // auto v = u->to[i]->node;
            // if (v->numFrom() != 1) {
            //     for (int j = 0; j < 4; ++j)
            //         if (v->from[j] != nullptr && v->from[j]->node == u) {
            //             v->from[j] = nullptr;
            //         }
            //     break;
            // }
            // u = v;
        }
        // printf("Walk through %d cvg = %d len = %d:", e->node->id, e->cvg, (int)path.size());
        // printf("%s", string(path.begin(), path.end()).c_str());
        // for (auto&& t : path)
        //     printf(" %d", t);
        // printf("\n");
    }
    void removeBubble()
    {
        int numBubble = 0;
        for (auto&& u : nodes) {
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    auto x = u->to[i];
                    auto y = u->to[j];
                    if (i != j && x != nullptr && y != nullptr) {
                        if ((double)x->cvg / y->cvg < 0.7) {
                            vector<pair<DBGNode*, int>> xpath, ypath;
                            xpath.push_back(make_pair(u, i));
                            ypath.push_back(make_pair(u, j));
                            walkThroughBubble(x->node, xpath);
                            walkThroughBubble(y->node, ypath);

                            auto tgt_x = xpath.back().first->to[xpath.back().second]->node;
                            auto tgt_y = ypath.back().first->to[ypath.back().second]->node;

                            if (tgt_x != tgt_y || xpath.size() < cfg.k / 2) {
                                printf("Walk x, target at %d len = %d cvg = %d\n", tgt_x->id, (int)xpath.size(), x->cvg);
                                for (auto&& t : xpath) {
                                    printf("%c", i2c(t.second));
                                }
                                puts("");
                                printf("Walk y, target at %d len = %d cvg = %d\n", tgt_y->id, (int)ypath.size(), y->cvg);
                                for (auto&& t : ypath) {
                                    printf("%c", i2c(t.second));
                                }
                                puts("");
                                continue;
                            }

                            ++numBubble;

                            auto lst = xpath.back().first;
                            auto tgt = lst->to[xpath.back().second]->node;
                            for (int k = 0; k < 4; ++k) {
                                if (tgt->from[k] != nullptr && tgt->from[k]->node == lst) {
                                    tgt->from[k] = nullptr;
                                }
                            }
                            u->to[i]->node->removed = true;
                            u->to[i] = nullptr;
                        }
                    }
                }
            }
        }
        printf("%d bubbles detected\n", numBubble);
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
            if (to != nullptr && !to->removed) {
                maxSucc = max(maxSucc, dfsExtLen(to->node->id));
            }
        }
        extLen[currId] += maxSucc;
        return extLen[currId];
    }
    DBGNode* calcExtLen()
    {
        extLen.resize(nodes.size());
        for (auto&& i : extLen)
            i = 0;
        int maxLen = 0;
        DBGNode* u;
        for (auto i = 0u; i < nodes.size(); ++i) {
            int len = dfsExtLen(i);
            if (maxLen < len) {
                maxLen = len;
                u = nodes[i];
            }
        }
        printf("Extend length: %d\n", maxLen);
        fflush(stdout);
        return u;
    }
    void followMaxExtLen(DBGNode* curr, vector<char>& path)
    {
        if (extLen[curr->id] != 1) {
            for (int i = 0; i < 4; ++i) {
                if (curr->to[i] != nullptr && !curr->to[i]->removed) {
                    auto next = curr->to[i]->node;
                    if (extLen[next->id] == extLen[curr->id] - 1) {
                        path.push_back(i2c(i));
                        curr->to[i]->removed = true;
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
    
    vector<Genome> exportPathsLengthFirst()
    {
        vector<Genome> res;
        while (true) {
            auto u = calcExtLen();
            if (extLen[u->id] > cfg.minOutputLength) {
                auto path = u->getRawData();
                followMaxExtLen(u, path);
                res.push_back(string(path.begin(), path.end()));
            } else {
                break;
            }
        }
        // for (auto&& u : nodes) {
        //     if (u->numFrom() == 0) {
        //         auto path = u->getRawData();
        //         auto calcExtLen();
        //         followMaxExtLen(u, path);
        //         res.push_back(string(path.begin(), path.end()));
        //     }
        // }
        return res;
    }
    void followLoopFirst(DBGNode* u, vector<char>& path)
    {
        bool hasLoop = false;
        for (int i = 0; i < 4; ++i) {
            auto v = u->to[i];
            if (v != nullptr && v->leadToLoop != 0 && !v->visited) {
                path.push_back(i2c(i));
                v->visited = true;
                followLoopFirst(v->node, path);
                hasLoop = true;
            }
        }
        if (!hasLoop) {
            for (int i = 0; i < 4; ++i) {
                auto v = u->to[i];
                if (v != nullptr && v->leadToLoop == 0) {
                    path.push_back(i2c(i));
                    followLoopFirst(v->node, path);
                    break;
                }
            }
        }
    }
    vector<Genome> exportPaths()
    {
        findLoop();
        vector<Genome> res;
        for (auto&& u : nodes) {
            if (u->numFrom() == 0) {
                auto path = u->getRawData();
                followLoopFirst(u, path);
                res.push_back(string(path.begin(), path.end()));
            }
        }
        return res;
    }
};
