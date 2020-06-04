#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Random.h"
#include "constants.h"
#include "kernel.h"

static void print_banner() {
  printf(
      "**                                                              **\n");
  printf(
      "** SciMark4 Numeric Benchmark, see http://math.nist.gov/scimark **\n");
  printf(
      "** for details. (Results can be submitted to pozo@nist.gov)     **\n");
  printf(
      "**                                                              **\n");
}

int main(int argc, char* argv[]) {
  /* default to the (small) cache-contained version */

  int itter = RESOLUTION_DEFAULT;
  Random R = new_Random_seed(RANDOM_SEED);

  if (argc > 1) {
    if (strcmp(argv[1], "-help") == 0 || strcmp(argv[1], "-h") == 0) {
      fprintf(stderr, "Usage: [itter]\n");
      exit(0);
    } else {
      itter = atoi(argv[1]);
    }
  }

  print_banner();
  /* print out results  */
  kernel_measure(itter, R);
  printf("\n");

  Random_delete(R);

  return 0;
}
