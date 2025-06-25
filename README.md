# Paper-Reproduce
In this repository, I will put the code and pdf files of my reproduced papers

## 0x00 
### d-MUL: Optimizing and Implementing a Multidimensional Scalar Multiplication Algorithm over Elliptic Curves
本文旨在探讨多维标量点乘算法d-MUL能否实现高效运算。已知d-MUL算法需要频繁进行高成本的矩阵运算和内存访问。论文第一部分推导了关于d-MUL加法链结构与构建的若干理论成果，这些理论发现本身就具有重要价值。第二部分则运用这些理论成果，提出了d-MUL的优化变体。我们的实现结果表明，d-MUL在维度d较小时具有很好的实用性，作为一种值得深入探索的算法，它在并行实现和密码学应用领域仍具有广阔的研究前景。
---
## 0x01
### A New Encoding Algorithm for a Multidimensional Version of the Montgomery Ladder
我们提出了一种新的编码算法，用于同时进行差分多维标量点乘算法 d-MUL。已有的编码算法在其高效和安全的实现方面存在重大缺陷。在2018年的一篇论文中，这些缺陷中的一些得到了避免，但代价是失去了点乘算法的一般功能。在本文中，我们解决了这些问题。我们的新编码算法以标量的二进制表示作为输入，构建一个紧凑的二进制序列和一个排列，这明确地确定了要在 d-MUL 中执行的组操作的规律序列。我们的算法简单地在标量上滑动大小为二的窗口，且效率极高。因此，在保持 d-MUL 完整的一般性的同时，我们成功消除了原始编码算法中的递归整数矩阵计算。我们还预计，我们的新编码算法将使 d-MUL 的实现更加容易，以恒定时间完成。我们的结果可以被视为将一维 Montgomery 梯子的高效和全面推广到任意维度。
