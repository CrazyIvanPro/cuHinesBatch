# QA

Q: 请问一下project里对串并行的精度差有多大的要求呢

A: 没有确切要求，可以自己在报告里写明实现的误差情况

Q: 预处理之类的复杂度是否计入？

按照题目上的描述 运行时间可仅包含Hines算法时间 那么预处理是否允许？

A: 允许预处理，只要不计算hines算法中需要计算的值都可以。

Q: 那么比如建立树形结构这种，也不计算到时间当中吗?

A: 以及想请问一下效率具体是怎么比较的...我感觉改动越多前面n小的点反而会更慢?

Q: 只要不计算hines算法中需要计算的值的处理操作都可以不计入。

A: 性能比较会重点参考后面四个规模较大的case。

Q: 请问proj里面给出的case全部都是树吗? 

A: 请仔细看Description文件，全部是树结构。

Q: 大作业要把输出结果也上传吗?

A: 去掉data、presult和sresult目录, 即上传你的串行、并行代码与报告即可。

Q: 请问助教如果使用openmp的话有必要在服务器上跑出结果吗，用自己电脑上的测试结果可以吗

A: 也可以，但最后评测性能以服务器为准，如果方便的话可以用wpn试一下。

Q: 哇最后还要性能比赛嘛

A: 大作业里面会有一个 占比不大 的性能分，应该会根据大家并行版本达到的性能 划分几个档次来给分。

Q: 以及 如果想了多种并行方法 最后提交代码只能提交一种吗?

A: 在报告说一下，从中优选一种。在project的目录中可以写一个README说一下以哪个版本为准。

Q: 怎样验证串行代码是正确的

A: 算法本身比较简单，根据伪码复现即可。重点需要参考数据的读取顺序、下标的边界等。在教学网的大作业一栏中也同步更新了一个助教写的参考程序的可执行文件（串行代码占一定的分数），分为Linux和Windows，供有需要的同学参考。

Q: 大作业要把输出结果也上传吗

A: 去掉data、presult和sresult目录, 即上传你的串行、并行代码与报告即可。



 



