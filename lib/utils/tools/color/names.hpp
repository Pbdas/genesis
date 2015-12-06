#ifndef GENESIS_UTILS_TOOLS_COLOR_NAMES_H_
#define GENESIS_UTILS_TOOLS_COLOR_NAMES_H_

/**
 * @brief
 *
 * @file
 * @ingroup utils
 */

#include "utils/tools/color.hpp"

#include <array>
#include <utility>
#include <string>

namespace genesis {

// =================================================================================================
//     Color Names
// =================================================================================================

struct ColorNames
{
    static const std::array<std::pair<std::string, Color>, 140> map;
};

} // namespace genesis

#endif // include guard
