# Searching for combinatorial designs based on Latin Squares

### Directories overview

/sat_enumeration_brown_dls - enumerate horizontally symmetric row-inverse and vertically symmetric column-inverse diagonal Latin squares of a given order via an AllSAT solver.

/solve_sat_cms - for a given CNF that encodes the search for two orthogonal Latin squares, a set of cells mapping schemas (CMS), and a set of partial Latin squares, for the Cartesian product of CMSs and partial Latin squares, add to the CNF
the CMS constraints, substitute the partial Latin square, and enumerate all solutions via an AllSAT solver.

/solve_sat_cms_mpi - an MPI version of solve_sat_cms.

/enumerate_brown_dls - enumerate orthogonal mates for a horizontally symmetric row-inverse diagonal Latin square of order 10 from the paper 'On the Construction of Triples of Diagonal Latin Squares of Order 10', form all triples based on them, and result in a triple with the highest orthogonality characteristic.

/dls_main_class_enumeration - for a given order n, a file with all diagonal Latin squares of order n normalized by the main diagonal, and a file of all ESODLS CMS, generate all main classes and all ESODLS CMS of order n and compare them with that from the file.

/dlx_mols - for a given file with Latin squares, find all their orthogonal mates by the DLX algorithm.
