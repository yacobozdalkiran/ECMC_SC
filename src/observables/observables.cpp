//
// Created by ozdalkiran-l on 1/8/26.
//

#include "observables.h"

// Computation of mean plaquette with halos embedded in field (needs field halos exchange first)
double mean_plaquette(const GaugeField& field, const GeometryCB& geo) {
    double sum = 0.0;
    SU3 U1, U2, U3, U4;
    for (size_t site = 0; site < geo.V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            for (int nu = mu + 1; nu < 4; nu++) {
                U1 = field.view_link_const(site, mu);
                U2 = field.view_link_const(geo.get_neigh(site, mu, up), nu);
                U3 = field.view_link_const(geo.get_neigh(site, nu, up), mu).adjoint();
                U4 = field.view_link_const(site, nu).adjoint();
                sum += (U1 * U2 * U3 * U4).trace().real() / 3.0;
            }
        }
    }
    return sum/((double)geo.V*6);
}

// Computation of the spatially averaged Polyakov loop (temporal direction mu = 3)
// P = 1/L^3 sum_x 1/3 Tr prod_{t=0}^{L-1} U_3(x, t)
double mean_polyakov_loop(const GaugeField& field, const GeometryCB& geo) {
    Complex sum(0.0, 0.0);
    SU3 P;
    size_t V3 = static_cast<size_t>(geo.L) * geo.L * geo.L;
    for (size_t site0 = 0; site0 < V3; site0++) {  // sites of the t = 0 slice
        P = SU3::Identity();
        size_t site = site0;
        for (int t = 0; t < geo.L; t++) {
            P *= field.view_link_const(site, 3);
            site = geo.get_neigh(site, 3, up);
        }
        sum += P.trace() / 3.0;
    }
    return std::abs(sum / static_cast<double>(V3));
}
