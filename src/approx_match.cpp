#include <string>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <cstring>
using namespace std;

#ifndef __APPROX
#define __APPROX

// O(n^2) to compute the minimal edit distance between two strings
int editDist(string a, string b)
{
    int n = a.length();
    int m = b.length();
    int f[n + 1][m + 1];
    memset(f, 0x3f, sizeof f);
    f[0][0] = 0;
    for (int i = 1; i <= n; ++i) {
        f[i][0] = i;
    }
    for (int j = 1; j <= m; ++j) {
        f[0][j] = j;
    }
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

pair<int, int> matchBF(string&S, string&T)
{
    int n = S.length(), m = T.length();
    int minDis = 0x3f3f3f3f, pos;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j <= n; ++j) {
            int d = editDist(S.substr(i, j - i), T);
            if (d < minDis) {
                minDis = d;
                pos = j;
            }
        }
    }
    printf("%d %d\n", minDis, pos);
    return make_pair(minDis, pos);
}


pair<int, int> match(string&S, string&T)
{
    int n = S.length(), m = T.length();
    int e[m + 1][n + 1], f[m + 1][n + 1];
    memset(e, 0x3f, sizeof e);
    memset(f, -1, sizeof f);
    for (int i = 0; i <= n; ++i)
        e[0][i] = 0;
    for (int i = 0; i <= m; ++i)
        e[i][0] = i;
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (e[i - 1][j] + 1 < e[i][j]) {
                e[i][j] = e[i - 1][j] + 1;
                f[i][j] = 0;
            }
            if (e[i][j - 1] + 1 < e[i][j]) {
                e[i][j] = e[i][j - 1] + 1;
                f[i][j] = 1;
            }
            if (T[i - 1] == S[j - 1]) {
                if (e[i - 1][j - 1] < e[i][j]) {
                    e[i][j] = e[i - 1][j - 1];
                    f[i][j] = 2;
                }
            } else {
                if (e[i - 1][j - 1] + 1 < e[i][j]) {
                    e[i][j] = e[i - 1][j - 1] + 1;
                    f[i][j] = 3;
                }
            }
        }
    }
    int minDis = 0x3f3f3f3f;
    int pos = -1;
    // for (int i = 0; i <= m; ++i) {
    //     for (int j = 0; j <= n; ++j) {
    //         printf("%d ", e[i][j]);
    //     }
    //     puts("");
    // }
    for (int i = 0; i <= n; ++i) {
        if (e[m][i] < minDis) {
            minDis = e[m][i];
            pos = i;
        }
    }
    printf("%d %d\n", minDis, pos);
    return make_pair(minDis, pos);
}



#endif