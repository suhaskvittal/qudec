/*
    author: Suhas Vittal
    date:   19 September 2025
*/

#include "argparse.h"

#include <iostream>

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

void
_die_with_error(std::string msg, std::string usage)
{
    std::cerr << usage << "\n";
    throw std::runtime_error(msg);
}

void
ARGPARSE::parse(int argc, char** argv)
{
    // define the usage string:
    std::string usage = "usage: " + std::string{argv[0]} + " " + usage_strm.str() + " [options]"
                        + "\n\nOPTIONS ---------------------------------------\n"
                        + options_strm.str();

    // now read the command line arguments:
    size_t required_idx{0};
    for (size_t i = 1; i < argc; ++i)
    {
        std::string x{argv[i]};
        
        // first check if we match a help flag/option:
        if (x == "-h" || x == "--help")
        {
            std::cout << usage << "\n";
            exit(0);
        }
    
        if (required_idx < required_arguments.size())
        {
            // now, make sure that this argument is not an option:
            const auto& [name, description, ptr, type] = required_arguments[required_idx];
            if (x.front() == '-')
            {
                _die_with_error("expected required argument `" 
                                        + std::string{name} + "` but got option `" + x + "`", usage);
            }
        
            read_argument_and_write_to_ptr(x, ptr, type);
            required_idx++;
        }
        else
        {
            if (x.front() != '-')
                _die_with_error("expected optional argument but got `" + x + "`", usage);
            
            bool is_option = (x[1] == '-');

            auto opt_it = std::find_if(optional_arguments.begin(), optional_arguments.end(),
                                        [&x, is_option] (const auto& arg) 
                                        { 
                                            return is_option ? arg.full_name == x : arg.flag_name == x; 
                                        });
            if (opt_it == optional_arguments.end())
                _die_with_error("unknown optional argument: " + x, usage);

            const auto& [flag_name, full_name, description, ptr, type] = *opt_it;
            if (type == ARGPARSE::TYPE_INFO::FLAG)
                *static_cast<bool*>(ptr) = true;
            else
                read_argument_and_write_to_ptr(std::string{argv[++i]}, ptr, type);
        }
    }

    if (required_idx < required_arguments.size())
    {
        _die_with_error("expected "
                + std::to_string(required_arguments.size() - required_idx) + " more required arguments",
                usage);
    }
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

template <class T> void
_read_helper(std::string arg, void* ptr)
{
    T* typed_ptr = static_cast<T*>(ptr);

    if constexpr (std::is_same<T, std::string>::value)
        *typed_ptr = arg;
    else if constexpr (std::is_same<T, int64_t>::value)
        *typed_ptr = std::stoll(arg);
    else if constexpr (std::is_same<T, double>::value)
        *typed_ptr = std::stod(arg);
    else
        throw std::runtime_error("argparse: unexpected type: " + std::string(typeid(T).name()));
}

void
read_argument_and_write_to_ptr(std::string arg, void* ptr, ARGPARSE::TYPE_INFO type)
{
    switch (type)
    {
    case ARGPARSE::TYPE_INFO::STRING:
        _read_helper<std::string>(arg, ptr);
        break;

    case ARGPARSE::TYPE_INFO::INT:
        _read_helper<int64_t>(arg, ptr);
        break;

    case ARGPARSE::TYPE_INFO::FLOAT:
        _read_helper<double>(arg, ptr);
        break;

    default:
        throw std::runtime_error("flag is unexpected -- should be resolved earlier.");
    }
}


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////