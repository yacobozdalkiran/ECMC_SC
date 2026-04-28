#include <iostream>

#include "../gauge/GaugeField.h"
#include "../observables/observables.h"
#include "../io/ildg.h"

void check(const GaugeField& field, const GeometryCB& geo) {
    double res = 0.0;
    for (size_t site = 0; site < geo.V; site++) {
        for (int mu = 0; mu < 4; mu++) {
            double temp =abs(field.view_link_const(site, mu).determinant()-1.0);
            if (temp>res) res = temp;
        }
    }
    std::cout << "Max |det(U)-1| = " << res <<"\n";
}

int main() {
    int L = 8;
    GeometryCB geo(L);
    GaugeField e5(geo);
    GaugeField e6(geo);
    GaugeField h5(geo);
    GaugeField h6(geo);
    read_ildg(e5, geo, "e5", "unitarity");
    read_ildg(e6, geo, "e6", "unitarity");
    read_ildg(h5, geo, "h5", "unitarity");
    read_ildg(h6, geo, "h6", "unitarity");
    check(e5, geo);
    check(e6, geo);
    check(h5, geo);
    check(h6, geo);
    std::cout << "P = " << mean_plaquette(e5, geo);
    std::cout << "P = " << mean_plaquette(e6, geo);
    std::cout << "P = " << mean_plaquette(h5, geo);
    std::cout << "P = " << mean_plaquette(h6, geo);
    return 0;
}
