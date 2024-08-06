/* Title: Getting Facets
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/

#include "utils.h"
#include "output.h"
#include "fractions.h"

scalar f[], f1[], f2[];
char filename[80];
bool includeCoat;

int main(int a, char const *arguments[]){
  sprintf(filename, "%s", arguments[1]);

  // Compare the string from command line argument to "true"
  if (strcmp(arguments[2], "true") == 0) {
    includeCoat = true;
    // fprintf(ferr, "includeCoat\n");
  } else {
    includeCoat = false;
    // fprintf(ferr, "not includeCoat\n");
  }

  restore (file = filename);

  if (includeCoat == true){
    foreach(){
      f[] = f1[];
    }
  } else {
    foreach(){
      f[] = f2[];
    }
  }

  FILE * fp = ferr;
  output_facets(f,fp);
  fflush (fp);
  fclose (fp);
}
