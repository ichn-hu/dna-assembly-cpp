#include <algorithm>
#include <chrono>
#include <cstdio>
#include <string>
#include <thread>
#include <vector>
using namespace std;

string hey[] = {
    "(     )",
    "(･    )",
    "(ω･   )",
    "(･ω･  )",
    "(´･ω･`)",
    "( ´･ω･)",
    "(  ´･ω)",
    "(   ´･)",
    "(    ´)",
};

int main()
{
    int i = 0;
    int lastLen = 0;
    int rep = 20;
    while (true) {
        for (int k = 0; k < rep; ++k)
            for (int j = 0; j < lastLen; ++j) {
                printf("\b");
            }
        for (int k = 0; k < rep; ++k)
            printf("%s ", hey[i].c_str());
        fflush(stdout);
        
        lastLen = (hey[i].length() + 1) * rep;
        i = (i + 1) % 9;
        std::chrono::milliseconds timespan(75);
        std::this_thread::sleep_for(timespan);
    }
    return 0;
}
