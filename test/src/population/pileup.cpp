/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2020

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

#include "src/common.hpp"

#include "genesis/population/formats/simple_pileup_reader.hpp"
#include "genesis/utils/text/string.hpp"

using namespace genesis::population;
using namespace genesis::utils;

TEST( Pileup, SimpleReader )
{
    // Skip test if no data availabe.
    NEEDS_TEST_DATA;
    std::string const infile = environment->data_dir + "population/example.pileup";

    auto reader = SimplePileupReader();
    auto records = reader.read( from_file( infile ));

    std::vector<char> ref_bases = { 'T', 'T', 'T', 'A', 'G', 'T', 'G', 'C' };

    ASSERT_EQ( 8, records.size() );
    for( size_t i = 0; i < records.size(); ++i ) {
        EXPECT_EQ( "seq1", records[i].chromosome );
        EXPECT_EQ( 272 + i, records[i].position );
        EXPECT_EQ( ref_bases[i], records[i].reference_base );

        ASSERT_EQ( 1, records[i].samples.size() );

        // LOG_DBG << i;
        // LOG_DBG1 << records[i].samples[0].read_bases;
        // LOG_DBG1 << join( records[i].samples[0].phred_scores );
    }

    EXPECT_EQ( "tTTTTTTttTtTtTTTtttTtTTT", records[0].samples[0].read_bases );
    EXPECT_EQ( "NNTTTTttTtTtTTTtttTtTTA",  records[1].samples[0].read_bases );
    EXPECT_EQ( "tTTT**ttTtTtTTTtttTtTTT",  records[2].samples[0].read_bases );
    EXPECT_EQ( "aAAAAaaAaAaAAAaaaAaAAAA",  records[3].samples[0].read_bases );
    EXPECT_EQ( "GGGTggGgGgGGGgggGgGGGG",   records[4].samples[0].read_bases );
    EXPECT_EQ( "TTTTttTtTtTCTtttTtTTGT",   records[5].samples[0].read_bases );
    EXPECT_EQ( "GGGGggGgGgGGGgggGgGGGGG",  records[6].samples[0].read_bases );
    EXPECT_EQ( "ACCTccCcC<><>cccCcCCCCC",  records[7].samples[0].read_bases );

    EXPECT_EQ(
        std::vector<unsigned char>({
            27, 27, 27, 10, 26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 28, 27, 26, 27, 26, 22, 27, 5
        }),
        records[0].samples[0].phred_scores
    );
    EXPECT_EQ(
        std::vector<unsigned char>({
            27, 27, 27, 26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 18, 27, 28, 27, 27, 27, 26, 27, 27, 10
        }),
        records[1].samples[0].phred_scores
    );
    EXPECT_EQ(
        std::vector<unsigned char>({
            22, 27, 22, 26, 27, 26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 28, 27, 26, 27, 26, 27, 27, 21
        }),
        records[2].samples[0].phred_scores
    );
    EXPECT_EQ(
        std::vector<unsigned char>({
            27, 10, 26, 24, 9, 27, 27, 27, 27, 27, 27, 27, 27, 27, 28, 27, 27, 25, 26, 27, 27, 27, 27
        }),
        records[3].samples[0].phred_scores
    );
    EXPECT_EQ(
        std::vector<unsigned char>({
            18, 18, 26, 10, 27, 27, 22, 28, 22, 27, 27, 22, 27, 5, 27, 27, 16, 26, 27, 27, 21, 27
        }),
        records[4].samples[0].phred_scores
    );
    EXPECT_EQ(
        std::vector<unsigned char>({
            10, 22, 27, 26, 27, 27, 27, 27, 27, 27, 27, 5, 27, 28, 27, 27, 25, 26, 27, 27, 5, 27
        }),
        records[5].samples[0].phred_scores
    );
    EXPECT_EQ(
        std::vector<unsigned char>({
            4, 18, 23, 9, 27, 27, 26, 27, 22, 27, 27, 22, 27, 28, 27, 27, 27, 26, 27, 27, 27, 27, 27
        }),
        records[6].samples[0].phred_scores
    );
    EXPECT_EQ(
        std::vector<unsigned char>({
            26, 22, 20, 5, 27, 27, 27, 27, 27, 27, 27, 27, 27, 28, 27, 27, 27, 24, 27, 27, 25, 27, 27
        }),
        records[7].samples[0].phred_scores
    );

    EXPECT_EQ(  0, records[0].samples[0].a_count );
    EXPECT_EQ(  0, records[0].samples[0].c_count );
    EXPECT_EQ(  0, records[0].samples[0].g_count );
    EXPECT_EQ( 24, records[0].samples[0].t_count );
    EXPECT_EQ(  0, records[0].samples[0].n_count );
    EXPECT_EQ(  0, records[0].samples[0].d_count );
    EXPECT_EQ( 24, records[0].samples[0].read_coverage );
    EXPECT_EQ( 24, records[0].samples[0].nucleotide_count );

    EXPECT_EQ(  1, records[1].samples[0].a_count );
    EXPECT_EQ(  0, records[1].samples[0].c_count );
    EXPECT_EQ(  0, records[1].samples[0].g_count );
    EXPECT_EQ( 20, records[1].samples[0].t_count );
    EXPECT_EQ(  2, records[1].samples[0].n_count );
    EXPECT_EQ(  0, records[1].samples[0].d_count );
    EXPECT_EQ( 23, records[1].samples[0].read_coverage );
    EXPECT_EQ( 21, records[1].samples[0].nucleotide_count );

    EXPECT_EQ(  0, records[2].samples[0].a_count );
    EXPECT_EQ(  0, records[2].samples[0].c_count );
    EXPECT_EQ(  0, records[2].samples[0].g_count );
    EXPECT_EQ( 21, records[2].samples[0].t_count );
    EXPECT_EQ(  0, records[2].samples[0].n_count );
    EXPECT_EQ(  2, records[2].samples[0].d_count );
    EXPECT_EQ( 23, records[2].samples[0].read_coverage );
    EXPECT_EQ( 21, records[2].samples[0].nucleotide_count );

    EXPECT_EQ( 23, records[3].samples[0].a_count );
    EXPECT_EQ(  0, records[3].samples[0].c_count );
    EXPECT_EQ(  0, records[3].samples[0].g_count );
    EXPECT_EQ(  0, records[3].samples[0].t_count );
    EXPECT_EQ(  0, records[3].samples[0].n_count );
    EXPECT_EQ(  0, records[3].samples[0].d_count );
    EXPECT_EQ( 23, records[3].samples[0].read_coverage );
    EXPECT_EQ( 23, records[3].samples[0].nucleotide_count );

    EXPECT_EQ(  0, records[4].samples[0].a_count );
    EXPECT_EQ(  0, records[4].samples[0].c_count );
    EXPECT_EQ( 21, records[4].samples[0].g_count );
    EXPECT_EQ(  1, records[4].samples[0].t_count );
    EXPECT_EQ(  0, records[4].samples[0].n_count );
    EXPECT_EQ(  0, records[4].samples[0].d_count );
    EXPECT_EQ( 22, records[4].samples[0].read_coverage );
    EXPECT_EQ( 22, records[4].samples[0].nucleotide_count );

    EXPECT_EQ(  0, records[5].samples[0].a_count );
    EXPECT_EQ(  1, records[5].samples[0].c_count );
    EXPECT_EQ(  1, records[5].samples[0].g_count );
    EXPECT_EQ( 20, records[5].samples[0].t_count );
    EXPECT_EQ(  0, records[5].samples[0].n_count );
    EXPECT_EQ(  0, records[5].samples[0].d_count );
    EXPECT_EQ( 22, records[5].samples[0].read_coverage );
    EXPECT_EQ( 22, records[5].samples[0].nucleotide_count );

    EXPECT_EQ(  0, records[6].samples[0].a_count );
    EXPECT_EQ(  0, records[6].samples[0].c_count );
    EXPECT_EQ( 23, records[6].samples[0].g_count );
    EXPECT_EQ(  0, records[6].samples[0].t_count );
    EXPECT_EQ(  0, records[6].samples[0].n_count );
    EXPECT_EQ(  0, records[6].samples[0].d_count );
    EXPECT_EQ( 23, records[6].samples[0].read_coverage );
    EXPECT_EQ( 23, records[6].samples[0].nucleotide_count );

    EXPECT_EQ(  1, records[7].samples[0].a_count );
    EXPECT_EQ( 17, records[7].samples[0].c_count );
    EXPECT_EQ(  0, records[7].samples[0].g_count );
    EXPECT_EQ(  1, records[7].samples[0].t_count );
    EXPECT_EQ(  0, records[7].samples[0].n_count );
    EXPECT_EQ(  0, records[7].samples[0].d_count );
    EXPECT_EQ( 23, records[7].samples[0].read_coverage );
    EXPECT_EQ( 19, records[7].samples[0].nucleotide_count );
}
