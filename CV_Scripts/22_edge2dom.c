#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <igraph.h>
#include <ctype.h>
#include "uthash.h"
#include "mat.h"

#define MAXCHAR 1024

/* This code takes as input an edge list. It uses this to create a
   network and then outputs the IDs of the nodes that are in the
   largest component of the network. 

Author: Joshua L. Payne */
int main(int argc, char * argv[])
{  

  int N = atoi(argv[1]); /* number of edges */
  FILE *fin = fopen(argv[2],"r"); /* name of input file */
  FILE *fout1 = fopen(argv[3],"w"); /* name of output file */
  FILE *fout2 = fopen(argv[4],"w"); /* name of output file */

  int i;
  int id;
  int biggest;
  char line[MAXCHAR]; /* place holder for each line */
  
  igraph_t graph;
  igraph_vector_t edges;
  igraph_vector_ptr_t components;
  igraph_vector_t membership;
  igraph_vector_t csize;
  igraph_integer_t no;
  
  /* initialize the edge vector */
  igraph_vector_init(&edges, 2*N);
  
  /* add edges to the edge vector */
  i = 0;
  while( fgets(line, MAXCHAR, fin) != NULL ) {
    
    id = atoi(strtok(line, "\n"));
    
    igraph_vector_set(&edges, i, id);
    i++;

  }
  
  /* use the edge vector to create the graph */
  igraph_create(&graph, &edges, 0, 0);
  
  /* decompose the graph into its components */
  igraph_vector_ptr_init(&components, 0);
  igraph_decompose(&graph, &components, IGRAPH_WEAK, -1, 2);

  /* write the edge list of the dominant component. remember that the
  vertices get relabeled, so we need to also write these vertices'
  original ids */
  igraph_write_graph_edgelist(VECTOR(components)[0], fout1);

  /* determine which vertices are in which components */
  igraph_vector_init(&membership, 0);
  igraph_vector_init(&csize, 0);
  igraph_clusters(&graph, &membership, &csize, &no, IGRAPH_WEAK);

  biggest = 0;
  for( i = 0; i < no; i++ ) {

    if( VECTOR(csize)[i] > biggest ) {
      
      biggest = i;
      
    }
    
  }


  for(i = 0; i < igraph_vcount(&graph); i++) {

    if( VECTOR(membership)[i] == biggest ) {
    
      fprintf(fout2, "%d\n", i);

    }

  }
    
  fclose(fin);
  fclose(fout1);
  fclose(fout2);
  
  return 0;

}
