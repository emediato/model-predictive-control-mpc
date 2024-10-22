#include <stdio.h>
#define LIM 25
#define FALSE 0
#define TRUE 1
  
     
void print_verify_simmetric(int TAM, A[LIM][LIM]){
  
  int simetrica = 1;
  for (int i = 1; i < TAM && simetrica; i++){
      for (int j = 0; j < i; j++){
          if (A[i][j] != A[j][i]){
              simetrica = 0;
              break; // se já achou um diferente, sai do for
          }
      }
  }
  if (simetrica)
      printf("é simétrica\n");
  else printf("não é simétrica\n");

} 
