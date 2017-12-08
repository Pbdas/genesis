/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2017 Lucas Czech

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lucas.czech@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

/**
 * @brief
 *
 * @file
 * @ingroup sequence
 */

#include "genesis/sequence/functions/codes.hpp"

#include "genesis/utils/text/string.hpp"
#include "genesis/utils/tools/color.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <unordered_map>

namespace genesis {
namespace sequence {

// =================================================================================================
//     Name Lists
// =================================================================================================

static const std::unordered_map<char, std::string> nucleic_acid_code_to_name = {
    { 'A', "Adenine" },
    { 'C', "Cytosine" },
    { 'G', "Guanine" },
    { 'T', "Thymine" },
    { 'U', "Uracil" },

    { 'W', "Weak" },
    { 'S', "Strong" },
    { 'M', "aMino" },
    { 'K', "Keto" },
    { 'R', "puRine" },
    { 'Y', "pYrimidine" },

    { 'B', "not A" },
    { 'D', "not C" },
    { 'H', "not G" },
    { 'V', "not T" },

    { 'N', "any" },
    { 'O', "omitted" },
    { 'X', "masked" },
    { '.', "gap" },
    { '-', "gap" },
    { '?', "gap" }
};

static const std::unordered_map<char, std::string> amino_acid_code_to_name = {
    { 'A', "Alanine" },
    { 'B', "Aspartic acid or Asparagine" },
    { 'C', "Cysteine" },
    { 'D', "Aspartic acid" },
    { 'E', "Glutamic acid" },
    { 'F', "Phenylalanine" },
    { 'G', "Glycine" },
    { 'H', "Histidine" },
    { 'I', "Isoleucine" },
    { 'J', "Leucine or Isoleucine" },
    { 'K', "Lysine" },
    { 'L', "Leucine" },
    { 'M', "Methionine" },
    { 'N', "Asparagine" },
    { 'O', "Pyrrolysine" },
    { 'P', "Proline" },
    { 'Q', "Glutamine" },
    { 'R', "Arginine" },
    { 'S', "Serine" },
    { 'T', "Threonine" },
    { 'U', "Selenocysteine" },
    { 'V', "Valine" },
    { 'W', "Tryptophan" },
    { 'Y', "Tyrosine" },
    { 'Z', "Glutamic acid or Glutamine" },

    { 'X', "any" },
    { '*', "translation stop" },
    { '-', "gap" },
    { '?', "gap" }
};

// =================================================================================================
//     Color Lists
// =================================================================================================

static const std::map<char, std::string> nucleic_acid_text_colors_map = {
    { 'A', "Red" },
    { 'C', "Green" },
    { 'G', "Yellow" },
    { 'T', "Blue" },
    { 'U', "Blue" },

    { 'W', "DarkGray" },
    { 'S', "DarkGray" },
    { 'M', "DarkGray" },
    { 'K', "DarkGray" },
    { 'R', "DarkGray" },
    { 'Y', "DarkGray" },

    { 'B', "DarkGray" },
    { 'D', "DarkGray" },
    { 'H', "DarkGray" },
    { 'V', "DarkGray" },

    { 'N', "DarkGray" },
    { 'O', "DarkGray" },
    { 'X', "DarkGray" },
    { '.', "DarkGray" },
    { '-', "DarkGray" },
    { '?', "DarkGray" }
};

static const std::map<char, std::string> amino_acid_text_colors_map = {
    { 'A', "Blue" },
    { 'B', "DarkGray" },
    { 'C', "LightMagenta" },
    { 'D', "Magenta" },
    { 'E', "Magenta" },
    { 'F', "Blue" },
    { 'G', "LightRed" },
    { 'H', "Cyan" },
    { 'I', "Blue" },
    { 'J', "DarkGray" },
    { 'K', "Red" },
    { 'L', "Blue" },
    { 'M', "Blue" },
    { 'N', "Green" },
    { 'O', "DarkGray" },
    { 'P', "Yellow" },
    { 'Q', "Green" },
    { 'R', "Red" },
    { 'S', "Green" },
    { 'T', "Green" },
    { 'U', "DarkGray" },
    { 'V', "Blue" },
    { 'W', "Blue" },
    { 'Y', "Cyan" },
    { 'Z', "DarkGray" },

    { 'X', "DarkGray" },
    { '*', "DarkGray" },
    { '-', "DarkGray" },
    { '?', "DarkGray" }
};

static const std::map<char, utils::Color> nucleic_acid_colors_map = {
    { 'A', { 1.0, 0.0, 0.0 } },
    { 'C', { 0.0, 1.0, 0.0 } },
    { 'G', { 1.0, 1.0, 0.0 } },
    { 'T', { 0.0, 0.0, 1.0 } },
    { 'U', { 0.0, 0.0, 1.0 } },

    { 'W', { 0.376, 0.376, 0.376 } },
    { 'S', { 0.376, 0.376, 0.376 } },
    { 'M', { 0.376, 0.376, 0.376 } },
    { 'K', { 0.376, 0.376, 0.376 } },
    { 'R', { 0.376, 0.376, 0.376 } },
    { 'Y', { 0.376, 0.376, 0.376 } },

    { 'B', { 0.5, 0.5, 0.5 } },
    { 'D', { 0.5, 0.5, 0.5 } },
    { 'H', { 0.5, 0.5, 0.5 } },
    { 'V', { 0.5, 0.5, 0.5 } },

    { 'N', { 0.5, 0.5, 0.5 } },
    { 'O', { 0.5, 0.5, 0.5 } },
    { 'X', { 0.5, 0.5, 0.5 } },
    { '.', { 0.5, 0.5, 0.5 } },
    { '-', { 0.5, 0.5, 0.5 } },
    { '?', { 0.5, 0.5, 0.5 } }
};

static const std::map<char, utils::Color> amino_acid_colors_map = {
    { 'A', { 0.098, 0.500, 1.000 } },
    { 'B', { 0.376, 0.376, 0.376 } },
    { 'C', { 0.902, 0.500, 0.500 } },
    { 'D', { 0.800, 0.302, 0.800 } },
    { 'E', { 0.800, 0.302, 0.800 } },
    { 'F', { 0.098, 0.500, 1.000 } },
    { 'G', { 0.902, 0.600, 0.302 } },
    { 'H', { 0.098, 0.702, 0.702 } },
    { 'I', { 0.098, 0.500, 1.000 } },
    { 'J', { 0.376, 0.376, 0.376 } },
    { 'K', { 0.902, 0.200, 0.098 } },
    { 'L', { 0.098, 0.500, 1.000 } },
    { 'M', { 0.098, 0.500, 1.000 } },
    { 'N', { 0.098, 0.800, 0.098 } },
    { 'O', { 0.376, 0.376, 0.376 } },
    { 'P', { 0.800, 0.800, 0.000 } },
    { 'Q', { 0.098, 0.800, 0.098 } },
    { 'R', { 0.902, 0.200, 0.098 } },
    { 'S', { 0.098, 0.800, 0.098 } },
    { 'T', { 0.098, 0.800, 0.098 } },
    { 'U', { 0.376, 0.376, 0.376 } },
    { 'V', { 0.098, 0.500, 1.000 } },
    { 'W', { 0.098, 0.500, 1.000 } },
    { 'Y', { 0.098, 0.702, 0.702 } },
    { 'Z', { 0.376, 0.376, 0.376 } },

    { 'X', { 0.5, 0.5, 0.5 } },
    { '*', { 0.5, 0.5, 0.5 } },
    { '-', { 0.5, 0.5, 0.5 } },
    { '?', { 0.5, 0.5, 0.5 } }
};

// =================================================================================================
//     Ambiguity Lists
// =================================================================================================

static const std::unordered_map<char, std::string> nucleic_acid_ambiguity_char_map = {
    { 'A', "A" },
    { 'C', "C" },
    { 'G', "G" },
    { 'T', "T" },
    { 'U', "T" },

    { 'W', "AT" },
    { 'S', "CG" },
    { 'M', "AC" },
    { 'K', "GT" },
    { 'R', "AG" },
    { 'Y', "CT" },

    { 'B', "CGT" },
    { 'D', "AGT" },
    { 'H', "ACT" },
    { 'V', "ACG" },

    { 'N', "ACGT" },
    { 'O', "-" },
    { 'X', "-" },
    { '.', "-" },
    { '-', "-" },
    { '?', "-" }
};

static const std::unordered_map< std::string, char > nucleic_acid_ambiguity_string_map = {
    { "A", 'A' },
    { "C", 'C' },
    { "G", 'G' },
    { "T", 'T' },

    { "AT", 'W' },
    { "CG", 'S' },
    { "AC", 'M' },
    { "GT", 'K' },
    { "AG", 'R' },
    { "CT", 'Y' },

    { "CGT", 'B' },
    { "AGT", 'D' },
    { "ACT", 'H' },
    { "ACG", 'V' },

    { "ACGT", 'N' },
    { "-", '-' }
};

// =================================================================================================
//     Codes
// =================================================================================================

// ---------------------------------------------------------------------
//     Nucleic Acids
// ---------------------------------------------------------------------

std::string nucleic_acid_codes_plain()
{
    return "ACGT";
}

std::string nucleic_acid_codes_degenerated()
{
    return "WSMKRYBDHV";
}

std::string nucleic_acid_codes_undetermined()
{
    return "NOX.-?";
}

std::string nucleic_acid_codes_all()
{
    return nucleic_acid_codes_plain()
         + nucleic_acid_codes_degenerated()
         + nucleic_acid_codes_undetermined();
}

// ---------------------------------------------------------------------
//     Amino Acids
// ---------------------------------------------------------------------

std::string amino_acid_codes_plain()
{
    return "ACDEFGHIKLMNOPQRSTUVWY";
}

std::string amino_acid_codes_degenerated()
{
    return "BJZ";
}

std::string amino_acid_codes_undetermined()
{
    return "X*-?";
}

std::string amino_acid_codes_all()
{
    return amino_acid_codes_plain()
         + amino_acid_codes_degenerated()
         + amino_acid_codes_undetermined();
}

// ---------------------------------------------------------------------
//     Misc
// ---------------------------------------------------------------------

std::string normalize_codes( std::string const& alphabet )
{
    // Uppercase, sort, uniq the alphabet.
    auto normalized = utils::to_upper_ascii( alphabet );
    std::sort( normalized.begin(), normalized.end() );
    normalized.erase( std::unique( normalized.begin(), normalized.end() ), normalized.end() );
    return normalized;
}

// =================================================================================================
//     Color Codes
// =================================================================================================

std::map<char, std::string> nucleic_acid_text_colors()
{
    return nucleic_acid_text_colors_map;
}

std::map<char, std::string> amino_acid_text_colors()
{
    return amino_acid_text_colors_map;
}

std::map<char, utils::Color> nucleic_acid_colors()
{
    return nucleic_acid_colors_map;
}

std::map<char, utils::Color> amino_acid_colors()
{
    return amino_acid_colors_map;
}

// =================================================================================================
//     Translate Codes
// =================================================================================================

std::string nucleic_acid_name( char code )
{
    auto ucode = toupper(code);
    if( nucleic_acid_code_to_name.count( ucode ) == 0 ) {
        throw std::out_of_range( "Invalid nucleic acid code '" + std::string( 1, code ) + "'." );
    }
    return nucleic_acid_code_to_name.at( ucode );
}

std::string amino_acid_name( char code )
{
    auto ucode = toupper(code);
    if( amino_acid_code_to_name.count( ucode ) == 0 ) {
        throw std::out_of_range( "Invalid amino acid code '" + std::string( 1, code ) + "'." );
    }
    return amino_acid_code_to_name.at( ucode );
}

std::string nucleic_acid_ambiguities( char code )
{
    auto ucode = toupper(code);
    if( nucleic_acid_code_to_name.count( ucode ) == 0 ) {
        throw std::out_of_range( "Invalid nucleic acid code '" + std::string( 1, code ) + "'." );
    }
    return nucleic_acid_ambiguity_char_map.at( ucode );
}

char nucleic_acid_ambiguity_code( std::string codes )
{
    // Uppercase, sort, uniq the codes.
    auto tmp = utils::to_upper_ascii( codes );
    std::sort( tmp.begin(), tmp.end() );
    tmp.erase( std::unique( tmp.begin(), tmp.end() ), tmp.end() );

    if( nucleic_acid_ambiguity_string_map.count( tmp ) == 0 ) {
        throw std::out_of_range( "Invalid nucleic acid codes '" + codes + "'." );
    }
    return nucleic_acid_ambiguity_string_map.at( tmp );
}

} // namespace sequence
} // namespace genesis
