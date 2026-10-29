/*
    author: Suhas Vittal
    date:   5 September 2025
*/

#ifndef ARGPARSE_h
#define ARGPARSE_h

#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>
#include <sstream>
#include <type_traits>
#include <vector>

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

/*
    For simplicity, we implement using a builder design pattern:
*/
class ARGPARSE
{
public:
    enum class TYPE_INFO { STRING, INT, FLOAT, FLAG };

    struct required_argument_type
    {
        std::string_view name;
        std::string_view description;
        void* ptr;
        TYPE_INFO type;
    };

    struct optional_argument_type
    {
        std::string_view flag_name{""};  // i.e., '-v'
        std::string_view full_name{""};  // i.e., '--verbose'
        std::string_view description;
        void* ptr;
        TYPE_INFO type;
    };
private:
    std::vector<required_argument_type> required_arguments;
    std::vector<optional_argument_type> optional_arguments;

    std::stringstream usage_strm;
    std::stringstream options_strm;
public:
    ARGPARSE() =default;

    template <class T> ARGPARSE& required(std::string_view name, std::string_view description, T& ref);
    template <class T, class DT> ARGPARSE& optional(std::string_view flag_name, std::string_view full_name,
                                            std::string_view description, T& ref, DT default_value);
    void parse(int argc, char** argv);
};

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

void read_argument_and_write_to_ptr(std::string, void*, ARGPARSE::TYPE_INFO);

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

template <class T> constexpr ARGPARSE::TYPE_INFO
argparse_get_type_info()
{
    if constexpr (std::is_same<T, std::string>::value)
        return ARGPARSE::TYPE_INFO::STRING;
    else if constexpr (std::is_same<T, int64_t>::value)
        return ARGPARSE::TYPE_INFO::INT;
    else if constexpr (std::is_same<T, double>::value)
        return ARGPARSE::TYPE_INFO::FLOAT;
    else
        return ARGPARSE::TYPE_INFO::FLAG;
}

template <class T> constexpr void
argparse_check_valid_type()
{
    constexpr bool type_is_ok = std::is_same<T, std::string>::value
                                || std::is_same<T, int64_t>::value
                                || std::is_same<T, double>::value
                                || std::is_same<T, bool>::value;
    static_assert(type_is_ok,
        "invalid type for argparse, only valid types are std::string, int64_t, double, and bool");
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

template <class T> ARGPARSE&
ARGPARSE::required(std::string_view name, std::string_view description, T& ref)
{
    if (optional_arguments.size() > 0)
        throw std::runtime_error("required arguments must be added before optional arguments");

    argparse_check_valid_type<T>();

    required_arguments.push_back({name, description, static_cast<void*>(&ref),
                                    argparse_get_type_info<T>()});

    usage_strm << " <" << name << ">";
    options_strm << std::setw(72) << std::left << name
                << std::setw(80) << std::left << description
                << std::setw(8) << std::left << "string"
                << std::setw(24) << std::left << "required" << "\n";

    return *this;
}

template <class T, class DT> ARGPARSE&
ARGPARSE::optional(std::string_view flag_name,
                        std::string_view full_name,
                        std::string_view description,
                        T& ref,
                        DT default_value)
{
    argparse_check_valid_type<T>();

    ref = static_cast<T>(default_value);
    optional_arguments.push_back({flag_name, full_name, description, static_cast<void*>(&ref),
                                    argparse_get_type_info<T>()});

    std::string name_string;
    if (flag_name.empty())
        name_string = std::string{full_name};
    else if (full_name.empty())
        name_string = std::string{flag_name};
    else
        name_string = std::string{flag_name} + ", " + std::string{full_name};

    std::string type_string;
    if constexpr (std::is_same<T, std::string>::value)
        type_string = "string";
    else if constexpr (std::is_same<T, int64_t>::value)
        type_string = "int";
    else if constexpr (std::is_same<T, double>::value)
        type_string = "float";
    else if constexpr (std::is_same<T, bool>::value)
        type_string = "bool";

    options_strm << std::setw(72) << std::left << name_string
                << std::setw(80) << std::left << description
                << std::setw(8) << std::left << type_string
                << std::setw(24) << std::left << "optional, default: " << default_value
                << "\n";

    return *this;
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

namespace
{

inline void
_die_with_error(std::string msg, std::string usage)
{
    std::cerr << usage << "\n";
    throw std::runtime_error(msg);
}

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

} // anon

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

inline void
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

inline void
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

#endif  // ARGPARSE_h
