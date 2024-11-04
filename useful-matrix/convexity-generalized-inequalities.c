#include <stdio.h>
#define LIM 25
#define FALSE 0
#define TRUE 1

theta = 3.14;
/*
Convexity with respect to generalized inequalities

f : R^n -> R^m is K-convex if dominium f is convex and

f (theta * x (1 - theta)*y ) <_k theta * f(x) + (1-theta)*f(y) ;
for x, y \in dom f, 0 <= theta <= 1

ex f : S^m -> S^m, f(X) = X^2 is S^m + -convex

proof: for fixed z \in R^m, z'X^2z = norm(X_z) is convex in X, i.e.,

z'(theta*X (1 - theta)*Y)^2 * z <= theta* z' X^2z + (1 - theta)z'*Y^2*z

for X, Y \in S^m, 0 <= theta <= 1

therefore

(theta * X (1 - theta)*Y )^2 <= theta * x^2 + (1-theta)*Y^2 

*/

     
void verify_convexity_with_respect_to_generalized_inequalities(int TAM, A[LIM][LIM]){
  
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

/*

convex optimization problems are tractable

*/
