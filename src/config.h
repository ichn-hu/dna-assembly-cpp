#include <string>
using namespace std;
const string DATA_ROOT_PATH = "./data/";

struct Config {
    string root;
    string longPath;
    string shortPath1;
    string shortPath2;
    string resultPath;
    string patchesPath;
    // split reads into kmers
    int k;
    int minOutputLength;
    Config(string name, int k = 29, int minOutputLength = 5000)
        : k(k)
        , minOutputLength(minOutputLength)
    {
        root = DATA_ROOT_PATH;
        longPath = DATA_ROOT_PATH + name + "/long.fasta";
        shortPath1 = DATA_ROOT_PATH + name + "/short_1.fasta";
        shortPath2 = DATA_ROOT_PATH + name + "/short_2.fasta";
        resultPath = DATA_ROOT_PATH + name + "/result.fasta";
        patchesPath = DATA_ROOT_PATH + name + "/patches.txt";
    }
};

#ifdef DATA1
Config cfg("data1", 29);
#define CONFIG
#endif
#ifdef DATA2
Config cfg("data2", 31);
#define CONFIG
#endif
#ifdef DATA3
Config cfg("data3");
#define CONFIG
#endif
#ifdef DATA4
Config cfg("data4");
#define CONFIG
#endif

#ifndef CONFIG
Config cfg("data1");
#endif
