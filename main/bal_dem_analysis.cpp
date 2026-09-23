/*
 *  author: Suhas Vittal
 *  date:   21 September 2026
 * */

#include <decoder/bunchaluts/hypergraph.h>
#include <circuit_generator.h>

#include <argparse/argparse.h>
#include <stim.h>

#include <iostream>

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

int
main(int argc, char* argv[])
{
    int64_t d;
    double p;

    ARGPARSE()
        .optional("-d", "--distance", "Distance of studied code", d, 5)
        .optional("-p", "--physical-error-rate", "Physical error rate parameterizing SI1000 error model", p, 1e-3)
        .parse(argc,argv);

    // Build circuit
    auto circ = toric_si1000(d, 3*d, p, false);
    // convert circuit into DEM
    auto dem = stim::circuit_to_dem(circ, {.decompose_errors = true});
    
    // create hypergraph and analyze it
    auto G = decoder::bal::hg::from_toric_code_dem(dem);

    std::cout << "Detectors = " << G.N()
                << ", Edges = " << G.M()
                << "\n";
    
    std::cout << "\n==========================================\n\n";

    const size_t dpr = d*d;
    const size_t base = (3*d/2)*dpr;

    std::cout << "Base = " << base << ", detectors per round = " << dpr 
                << ", base degree = " << G.degree(base)
                << "\n";

    bool all_detectors_have_same_degree{true};
    for (size_t i = dpr; i < G.N()-dpr; i++)
        if (G.degree(i) != G.degree(base))
            all_detectors_have_same_degree = false;
    std::cout << "All detectors have same degree? (between " << dpr << " <= D < " << (G.N()-dpr) << ")?  --> "
                << all_detectors_have_same_degree << "\n";

    for (size_t k = 1; k <= d/2; k++)
    {
        auto cu = decoder::bal::hg::extract_error_cube(G, base, k);
        std::cout << k << "-cube:"
                << "\n\tN = " << cu.nodes.size()
                << "\n\tM = " << cu.edges.size()
                << "\n";
    }
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
