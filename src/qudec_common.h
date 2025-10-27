/*
 *  author: Suhas Vittal
 *  date:   12 October 2025
 */

#ifndef QUDEC_COMMON_h
#define QUDEC_COMMON_h

/*
 * Utility functions used in executables:
 * */

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

template <class T> void
print_stat(std::ostream& out, std::string name, T value)
{
    out << std::setw(64) << std::left << name;
    if constexpr (std::is_floating_point<T>::value)
    {
        if (value < 1e-3)
            out << std::setw(12) << std::right << std::scientific << std::setprecision(4) << value;
        else
            out << std::setw(12) << std::right << std::fixed << std::setprecision(8) << value;
    }
    else
    {
        out << std::setw(12) << std::right << value;
    }
    out << "\n";
}

template <class INTA, class INTB> double
fpdiv(INTA a, INTB b)
{
    return static_cast<double>(a) / static_cast<double>(b);
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

template <class IMPL> void
print_decoder_stats(std::ostream& out, const IMPL& dec)
{
}

template <class IMPL, class... IMPL_ARGS> DECODER_STATS
eval_decoder(const stim::Circuit& circuit, uint64_t num_trials, IMPL_ARGS... args)
{
    IMPL decoder(circuit, std::forward<IMPL_ARGS>(args)...);
    auto out = benchmark_decoder(circuit, decoder, num_trials);

    print_decoder_stats(std::cout, decoder);

    return out;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////

#endif  // QUDEC_COMMON_h
