/**
*  BP_algo_one_path.cpp
*  Used to solve the constrained longest path problem for an acyclic graph,
*  starting from given start nodes to given endnodes. Uses a cell array data structure for the
*  network.
*/

#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include <Math.h>
#include "matrix.h"
#include "mex.h"
#include <vector>
using namespace std;


/*
#include <cmath>
Gateway function

Inputs:
prhs[0] - a cell array where cell k is a vector corresponding to node k in the
graph and holds all the nodes that node k connects to.
prhs[1] - a cell array with the costs corresponding to the arcs described in prhs[0].
prhs[2] - total number of nodes in graph.
prhs[3] - a vector with the startnodes in a topological order (nodes from where paths can start in the acyclic network)
prhs[4] - a vector with the endnodes in a topological order (nodes from where paths can end in the acyclic network)
prhs[5] - maximum size of path
prhs[6] - index of the first endnode
prhs[7] - cell containing the parameters of the graph as returned by the function "GetEdgeletGraphParam.m"

Outputs:
plhs[0] - a vector with costs of all the best paths of size l <= maxlength
plhs[1] - list of predecessor nodes for the shortest path
*/


// Declaration of the function used by the Matlab mex function:
void path2TF(int *t, double *f, int *p, int J, int fmin, int fmax, int fsc, double delta_p, int length_p);

// Matlab mex function :
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray
	*prhs[])
{

	// -------------------------------------------------------------------
	// Initialization some variables :
	// -------------------------------------------------------------------

	double *arcs;         // the set of endpoints and costs of arcs emanating out of node k
	int nArcs;            // number of arcs
	double *arcCosts;     // the arc costs
	double *startNodes;   // starting nodes in the graph
	double *endNodes;     // end nodes in the graph
	int nStartNodes;      // length of the array start nodes
	int nEndNodes;        // length of the array end nodes 
	int nNodes;           // number of nodes in graph
	int maxLength;

	const mxArray *cellElementPtr; // Ptr on the arcs, or the corresponding costs, or on the graph parameters

	int startNode;
	int arcIndex;
	double tentativeDistance;
	int k;
	int l;

	int J;
	int fsc;
	int minfreq;
	int maxfreq;
	double delta_p;
	double C_alpha;
	double alpha;



	// -------------------------------------------------------------------
	// Getting the input data (communicating with Matlab) :
	// -------------------------------------------------------------------

	nNodes = (int)mxGetScalar(prhs[2]);
	startNodes = mxGetPr(prhs[3]); // We do not need the startNodes, we only need the number of startNodes
	nStartNodes = mxGetNumberOfElements(prhs[3]);
	endNodes = mxGetPr(prhs[4]);
	nEndNodes = mxGetNumberOfElements(prhs[4]);
	maxLength = (int)mxGetScalar(prhs[5]);
	int nLeftEndpoints;   // index of the first endnode
	nLeftEndpoints = (int)mxGetScalar(prhs[6]);
	
	
	// Getting the graph parameters :
	double *param;

	cellElementPtr = mxGetCell(prhs[7], 0);
	param          = mxGetPr(cellElementPtr);
	J              = (int)param[1];

	cellElementPtr = mxGetCell(prhs[7], 1);
	param          = mxGetPr(cellElementPtr);
	fsc            = (int)param[1];

	cellElementPtr = mxGetCell(prhs[7], 3);
	param          = mxGetPr(cellElementPtr);
	minfreq        = (int)param[0];
	maxfreq        = (int)param[1];

	cellElementPtr = mxGetCell(prhs[7], 4);
	param          = mxGetPr(cellElementPtr);
	delta_p        = (double)param[0];

	cellElementPtr = mxGetCell(prhs[7], 5);
	param          = mxGetPr(cellElementPtr);
	C_alpha        = param[0];
	alpha          = param[1];

	

	// ---------------------------------------------------------------------------------------------------------------------------
	// Initialization of matrix d  and pred containing respectively the distances of the nodes and the predecessors (with Matlab):
	// ---------------------------------------------------------------------------------------------------------------------------

	//// Allocation of memory for the matrix d and pred :
	int size_l = maxLength + 1;

	plhs[0] = mxCreateDoubleMatrix(size_l, nNodes, mxREAL);
	double *d;
	d = mxGetPr(plhs[0]);
	int size_d = size_l * nNodes;
	for (int i = 0; i<size_d; i++) {
		d[i] = 0;
	}

	plhs[1] = mxCreateDoubleMatrix(size_l, nNodes, mxREAL);
	double *pred;
	pred = mxGetPr(plhs[1]);
	for (int i = 0; i<size_d; i++) {
		pred[i] = 0;
	}

	int incr = 1;
	for (int i = 0; i < nStartNodes; i++){
		pred[i * size_l] = incr;
		incr++;
	}



	// --------------------------------------------------------------------
	// Allocation of the memory for t, f and path :
	// --------------------------------------------------------------------

	int* t    = (int *)malloc(3 * sizeof(int));
	double* f = (double *)malloc(3 * sizeof(int));
	int* path = (int *)malloc(3 * sizeof(int));

	double coeff1;
	double coeff2;
	double angle;
	double angle_max = 3.14 / 6;



	// -------------------------------------------------------------------
	// Start of the loop on the nodes of the graph (in topological order) :
	// -------------------------------------------------------------------

	for (startNode = 0; startNode < nLeftEndpoints; startNode++) {




		// -----------------------------------------------------------------------
		// Getting the arcs and its corresponding costs (given by the Matlab code)
		// -----------------------------------------------------------------------

		cellElementPtr = mxGetCell(prhs[0], startNode);

		if (cellElementPtr == NULL) { // the cell could be empty, so in that case we do nothing...
		continue;
		}

		// we reach this statement if the cell is not empty...

		// We get the arcs emanating from startNode :
		arcs = mxGetPr(cellElementPtr);
		nArcs = mxGetNumberOfElements(cellElementPtr);

		// We get the corresponding costs of these arcs :
		cellElementPtr = mxGetCell(prhs[1], startNode);
		arcCosts = mxGetPr(cellElementPtr);

	

		// -------------------------------------------------------------------
		// Start of the loop on the different arcs emanating from startNode :
		// -------------------------------------------------------------------
		for (k = 0; k < nArcs; k++) {



			// -------------------------------------------------------------------
			// Start of the loop on the different size of paths :
			// -------------------------------------------------------------------
			for (l = 1; l<maxLength + 1; l++) { // +1 to take care of zero length



				// -----------------------------------------------------------------------------------------------------
				// Restrictions on the best path depending on the size l of the path, and the angles between the nodes :
				// -----------------------------------------------------------------------------------------------------

				// We check out whether or not the path that we are computing makes sense:
				if (pred[l - 1 + startNode * size_l] == 0){ 
					continue; // Then the path does not make sens...
				}


				// We check out whether or not the angle between startNode and arcIndex is too large:
				arcIndex = (int)arcs[k] - 1;

				if (startNode > nStartNodes - 1) {
					path[0] = pred[l - 1 + startNode * size_l];
					path[1] = startNode;
					path[2] = arcIndex;
					path2TF(t, f, path, J, minfreq, maxfreq, fsc, delta_p, (int)3);
					coeff1 = (double)(f[1] - f[0]) / (t[1] - t[0]);
					coeff2 = (double)(f[2] - f[1]) / (t[2] - t[1]);
					angle  = (double)(coeff1 - coeff2) / (1 + coeff1 * coeff2);
					angle  = (double)atan(angle);

					if (angle > angle_max){
						continue; // Then the angle is too large...
					}

				}



				// -------------------------------------------------------------------
				// Updating the distance matrix d and pred if needed:
				// -------------------------------------------------------------------

				tentativeDistance = d[l - 1 + startNode * size_l] + arcCosts[k];

				if (fabs(d[l + arcIndex * size_l ]) <= fabs( tentativeDistance ) ) {
				d[l + arcIndex * size_l] = tentativeDistance; // change the value of the distance
				pred[l + arcIndex * size_l] = startNode + 1;  // change predecessor node
				}


			} // End of the loop on the different sizes of paths

		} // End of the loop on the different arcs emanating from startNode:

	} // End of the loop on the nodes of the graph (in a topological order)



	// -------------------------------------------------------------------
	// We display the distance and predecessor matrix (for testing the code) :
	// -------------------------------------------------------------------

	/*
	for (int i = 0; i < size_l; i++) // line i
	{
		for (int j = 0; j < nNodes; j++) // column j
		{
			printf(" %f", d[i][j]);
		}
		printf("\n");
	}

	printf("\n");
	printf("\n");

	for (int i = 0; i < size_l; i++) // line i
	{
		for (int j = 0; j < nNodes; j++) // column j
		{
			printf(" %d", pred[i][j]);
		}
		printf("\n");
	}
	*/

	// -------------------------------------------------------------------
	// Free memory:
	// -------------------------------------------------------------------

	free(t);
	free(f);
	free(path);


} // end of Matlab mex function







// -------------------------------------------------------------------
// Description of the function PATH2TF:
// -------------------------------------------------------------------

/*
Path2TF:
Returns the x and y indices
for the given edgelet path p as returned by
CalculateStatistic.

Inputs
p		array of the topological order of the nodes along a path
J       dyadic length of signal
fmin    width of the upper horizontal strip not computed by the algorithm, because too close to the upper edge (need to be an integer)
fmax    width of the bottom horizontal strip not computed by the algorithm, because too close to the bottom edge (need to be an integer)
fsc     the finest scale in the chirplet gragh (needs to be an integer: 0 <= csc <= fsc <= J )
Outputs
t       coordinates on the x-axis for the edgelet path
freq    coordinates on the y-axis for the edgelet path

See Also
- DisplayEdgelets, DisplayEdgeletPath
- CalculateStatistic
*/



void path2TF(int *t, double *f, int *p, int J, int fmin, int fmax, int fsc, double delta_p, int length_p)
{

	int nfreqs;
	nfreqs = (int)(fmax - fmin + 1) / delta_p;

	for (int i = 0; i < length_p; i++){

		f[i] = (p[i] - 1) % nfreqs; // modulo
		
		t[i] = (p[i] - f[i] - 1) / nfreqs;
		t[i] = t[i] * pow(2, (J - fsc));

		f[i] = f[i] * delta_p; // We need to change the value of the frequence if we are considering an sub - pixelic image

		f[i] = f[i] + fmin + 1 * delta_p;

	}

}