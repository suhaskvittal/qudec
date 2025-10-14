/*
    author: Suhas Vittal
    date:   5 September 2025
*/

#ifndef ARGPARSE_h
#define ARGPARSE_h

#include <algorithm>
#include <stdexcept>
#include <iomanip>
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

#endif  // ARGPARSE_h