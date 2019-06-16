#include "debrujin.cpp"

main()
{
    auto genomes = readFasta(cfg.shortPath1);
    auto tmp = readFasta(cfg.shortPath2);
    vector<Genome> tmp_rev_comp;
    for (auto&&t:tmp) tmp_rev_comp.push_back(t.reverse().compliment());
    genomes.insert(genomes.end(), tmp_rev_comp.begin(), tmp_rev_comp.end());
    // extendGenomes(genomes);
    cfg.minOutputLength = 400;
    auto graph = new DeBrujinGraph(genomes, cfg);
    graph->analysisNodeInOut();
    graph->removeBubble();
    // graph->findLoop();
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

    auto graph2 = new DeBrujinGraph(res, Config(cfg.root, 3000));
    graph2->analysis();
    graph2->removeBubble();
    res = graph2->exportPaths();
    writeFasta(res, cfg.resultPath, cfg.minOutputLength);
    // graph->output();
    return 0;
}
