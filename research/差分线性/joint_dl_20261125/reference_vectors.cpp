#include "../gaston_dl/include/gaston.hpp"
#include "../hsilu_dl/include/hsilu.hpp"
#include <iostream>
template<class S> void print(const S& s){for(auto x:s)std::cout<<x<<' ';}
int main(){
  for(int k=0;k<32;k++){
    auto x=gdl::sample_state(20261125,k);
    print(x);print(gdl::linear(x));print(hdl::linear(x));
    print(gdl::core(x,4,8,true));print(hdl::core(x,4));std::cout<<'\n';
  }
}
