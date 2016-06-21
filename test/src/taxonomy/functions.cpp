/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2016 Lucas Czech

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
 * @ingroup test
 */

#include "common.hpp"

#include "lib/taxonomy/formats/taxonomy_reader.hpp"
#include "lib/taxonomy/formats/taxscriptor_generator.hpp"
#include "lib/taxonomy/formats/taxscriptor_parser.hpp"
#include "lib/taxonomy/functions/taxonomy.hpp"
#include "lib/taxonomy/functions/taxscriptor.hpp"
#include "lib/taxonomy/taxon.hpp"
#include "lib/taxonomy/taxonomy.hpp"
#include "lib/taxonomy/taxscriptor.hpp"

using namespace genesis::taxonomy;

TEST( Taxonomy, Counts )
{
    // Skip test if no data availabe.
    NEEDS_TEST_DATA;
    std::string infile;

    auto reader = TaxonomyReader();
    reader.rank_field_position( 2 );

    // Read file.
    Taxonomy tax;
    infile = environment->data_dir + "taxonomy/tax_slv_ssu_123.1.unordered";
    EXPECT_NO_THROW( reader.from_file( infile, tax ));
    EXPECT_EQ( 32, total_taxa_count(tax) );
    EXPECT_TRUE( validate( tax ));

    // Level count.
    EXPECT_EQ( 1, taxa_count_at_level( tax, 0 ));
    EXPECT_EQ( 4, taxa_count_at_level( tax, 1 ));
    EXPECT_EQ( 5, taxa_count_at_level( tax, 2 ));
    EXPECT_EQ( 7, taxa_count_at_level( tax, 3 ));
    EXPECT_EQ( 6, taxa_count_at_level( tax, 4 ));
    EXPECT_EQ( 9, taxa_count_at_level( tax, 5 ));
    EXPECT_EQ( 0, taxa_count_at_level( tax, 6 ));

    // Rank count.
    EXPECT_EQ( 1, taxa_count_with_rank( tax, "Domain" ));
    EXPECT_EQ( 4, taxa_count_with_rank( tax, "Phylum" ));
    EXPECT_EQ( 5, taxa_count_with_rank( tax, "Class" ));
    EXPECT_EQ( 7, taxa_count_with_rank( tax, "Order" ));
    EXPECT_EQ( 6, taxa_count_with_rank( tax, "Family" ));
    EXPECT_EQ( 9, taxa_count_with_rank( tax, "Genus" ));
    EXPECT_EQ( 0, taxa_count_with_rank( tax, "Something" ));
}
