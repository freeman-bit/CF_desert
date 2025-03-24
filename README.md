# CF_desert_open
Here is a solver I selected from my private code repository, which I use to demonstrate my personal skills to potential employers.

This solver is used to solve scalar physical field interface discontinuities using a collocated grid finite volume scheme. 
The bidirectional physical fields employ a fully implicit solving format, and more important, the jump condition is applied implicitly which subjected to a tolerance with large time step. 
For example, the temperature field solving for a bi-directional interface with thermal resistance, the algorithm has second-order accuracy.
