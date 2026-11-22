#include "hsilu.h"
void hsilu_round(hsilu_state_t*, uint64_t);
extern "C" void verify_round(uint64_t* p, uint64_t c) { hsilu_state_t a; for(int i=0;i<5;i++)a.x[i]=p[i]; hsilu_round(&a,c); for(int i=0;i<5;i++)p[i]=a.x[i]; }
