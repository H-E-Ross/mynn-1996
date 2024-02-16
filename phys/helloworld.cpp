#include <algorithm>
#include <functional>
#include <stdio.h>

extern "C" void helloworld();

using std::bind;

void helloworld() {
   // printf() displays the string inside quotation
   printf("Hello, World!");
   return;
}
