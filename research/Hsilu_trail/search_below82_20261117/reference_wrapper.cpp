#include "../../hsilu_ref/hsilu.cpp"

extern "C" void reference_linear(const uint64_t* in,uint64_t* out){
    int inv[32];
    for(int x=0;x<32;x++){
        int y=0;for(int i=0;i<5;i++)y|=(((x>>i&1)^((1^(x>>((i+1)%5)&1))&(x>>((i+2)%5)&1)))<<i);
        inv[y]=x;
    }
    hsilu_state_t state={};
    for(int z=0;z<64;z++){
        int y=0;for(int i=0;i<5;i++)y|=(in[i]>>z&1)<<i;
        int x=inv[y];for(int i=0;i<5;i++)state.x[i]|=uint64_t(x>>i&1)<<z;
    }
    hsilu_round(&state,0);
    for(int i=0;i<5;i++)out[i]=state.x[i];
}
