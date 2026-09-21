/*
 *  author: Suhas Vittal
 *  date:   21 September 2026
 * */

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
    auto circ = toric_si1000(d, d, p, false);
    // convert circuit into DEM
    auto dem = stim::circuit_to_dem(circuit, {.decompose_errors = true});
    
    const size_t dpr = d*d;
    const size_t base = dpr;

    bal::Hypergraph
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
