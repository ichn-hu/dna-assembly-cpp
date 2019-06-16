#include "debrujin.cpp"

main()
{
    auto genomes = readFasta(cfg.shortPath1);
    auto tmp = readFasta(cfg.shortPath2);
    genomes.insert(genomes.end(), tmp.begin(), tmp.end());
    extendGenomes(genomes);
    auto graph = new DeBrujinGraph(genomes, cfg);
    graph->analysisNodeInOut();
    graph->removeBubble();
    // graph->analysis();
    // graph->removeBubble();
    // graph->analysis();
    // graph->removeBubble();

    // graph->removeBubble();
    // graph->analysisBranch();
    // graph->analysis();
    // auto res = graph->exportPaths();
    auto res = graph->exportPathsLengthFirst();
    extendGenomes(res);

    // auto graph2 = new DeBrujinGraph(res, Config(cfg.root, 3000));
    // graph2->analysis();
    // graph2->removeBubble();
    // res = graph2->exportPaths();
    writeFasta(res, cfg.resultPath, cfg.minOutputLength);
    // graph->output();
    return 0;
}
