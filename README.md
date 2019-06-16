---
typora-copy-images-to: ./imgs
---

[TOC]

# 宏基因组组装PJ

胡志峰 16307130177

## 妹她

mailto: wanjunpeng@foxmail.com

## 问题背景

### 评价指标

- Genome fraction(%):组装好的contig比对到参考参考基因组后占参考基因组总长度的比例。
- NGA50: 将组装好的contig比对按长度从大到小匹配到到参考基因组，当匹配总占比超过原基因组50%时该contig的长度
  对不同物种分别计算后取平均值。我明白了=。=
- Duplication ratio: 将组装好的contig比对回基因组后参考基因组每个位置的平均重复率（即可能有多个contig覆盖了同一段参考基因组的序列，该指标需要尽可能小）。
- Misassemblies: 组装好的contig中连续的片段映射到参考基因组中非连续的片段上,则判定该contig为misassembly(组装错误),最后统计所有misassembly的数目。
- Mismataches per 100kbp: 匹配回参考基因组后每100000个字符中与参考基因组不同的字符数（可以理解为拼接结果的错误率,参考测序的错误率描述）



## 参考资料

[^1]: De Bruijn graph <https://en.wikipedia.org/wiki/De_Bruijn_graph>
[^2]: De novo sequence assemblers <https://en.wikipedia.org/wiki/De_novo_sequence_assemblers>
[^3]: 字符串相似度之美（一）<https://zhuanlan.zhihu.com/p/20101194>
[^4]: <https://zhuanlan.zhihu.com/p/20102352>