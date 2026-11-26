// Columnwise HNW recurrence with complete-mask evaluation through L^T.
// The last chi is expanded using squared LAT entries before this readout.
#define main original_hnw_cli_main
#ifdef TARGET_GASTON
#include "../gaston_dl/src/hnw_endpoint.cpp"
#else
#include "../hsilu_dl/src/hnw_endpoint.cpp"
#endif
#undef main
int main(){
    HNW h(false);
    State d{},mask{};
    while(std::cin >> std::hex >> d[0]){
        for(int i=1;i<5;i++)std::cin >> d[i];
        for(int i=0;i<5;i++)std::cin >> mask[i];
        Vector q(64,Column(32));
        for(int j=0;j<64;j++){
            for(int v=0;v<32;v++)q[j].n[v]=parity(v&column(d,j))?-1:1;
            q[j]=h.nonlinear(q[j]);
        }
        std::vector<int> columns;
        for(int j=0;j<64;j++)if(column(mask,j))columns.push_back(j);
        Dyadic total;
        std::function<void(size_t,State,Dyadic)> visit;
        visit=[&](size_t k,State alpha,Dyadic weight){
            if(k==columns.size()){
#ifdef TARGET_GASTON
                const State beta=h.masks.transpose(alpha);
#else
                const State beta=transpose(alpha);
#endif
                for(int j=0;j<64;j++)weight=weight*q[j].get(column(beta,j));
                total+=weight;
                return;
            }
            const int j=columns[k],b=column(mask,j);
            for(int a=0;a<32;a++)if(h.tables.lat[b][a]){
                State next=alpha,part=embed(j,a);
                for(int i=0;i<5;i++)next[i]^=part[i];
                visit(k+1,next,weight*Dyadic(h.tables.lat[b][a]*h.tables.lat[b][a],10));
            }
        };
        visit(0,State{},Dyadic(1));
        std::cout<<std::dec;json(std::cout,total);std::cout<<'\n';
    }
    return 0;
}
