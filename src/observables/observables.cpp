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
