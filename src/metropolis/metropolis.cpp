#include "metropolis.h"


void metropolis::sweep(GaugeField& field, const GeometryCB& geo, double beta, const std::vector<SU3>& set, size_t& accepted, size_t& proposed, std::mt19937_64& rng){
    //Effectue un sweep metropolis avec n_hits hits à chaque lien.
    SU3 staple;
    SU3 Unew;
    accepted = 0;
    proposed = 0;
    int i_set = 0;
    int n_hits=1;
    std::uniform_int_distribution<int> index_set(0, static_cast<int>(set.size()) - 1);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    for (size_t site = 0; site < geo.V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            //On copie la staple et Uold
            field.compute_staple(geo, site, mu, staple);
            auto Umap = field.view_link(site, mu);
            SU3 Uold = Umap;

            for (int i =0; i < n_hits; i++) {
                //On choisit une matrice d'update dans set
                i_set = index_set(rng);

                //On copie Unew
                Unew = set[i_set] * Uold;

                //On propose le move
                double old_tr = (Uold * staple).trace().real();
                double new_tr = (Unew * staple).trace().real();
                double dS = -(beta/3.0) * (new_tr - old_tr);

                ++proposed;
                bool accept = false;
                if ((dS <= 0.0)||(unif(rng) < exp(-dS))) accept = true;

                if (accept) {
                    Umap.noalias() = Unew; //Unew et Umap ne se chevauchent pas -> safe
                    ++accepted;
                }
            }
        }
    }
}
