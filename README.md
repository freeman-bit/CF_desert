# CF_desert (public part)
Here is a solver I selected from my private code repository, which I use to demonstrate my personal skills to potential employers.

This solver is used to solve scalar physical field interface discontinuities using a collocated grid finite volume scheme baseon baasilisk, an open-source software with a quadtree/octree numerical structure.  
The bidirectional physical fields employ a fully implicit solving format, and more important, the sharp-jump condition is applied implicitly which subjected to a tolerance with large time step. 
For example, the temperature/concentration field solving for a bi-directional interface with thermal resistance, the algorithm has second-order accuracy.

If you would like to learn more details, please contact weibolin@stu.xjtu.edu.cn.
