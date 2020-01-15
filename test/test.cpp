#include "nuclide_class.h"
#include "configure.h"
#include "gtest/gtest.h"
#include "reactions.h"
#include "../extern/pugiData/pugixml.h"
#include "../extern/xtensor/include/xtensor/xview.hpp"
#include "../extern/xtensor/include/xtensor/xbuilder.hpp"
#include "../extern/xtensor/include/xtensor/xcomplex.hpp"

#include "../extern/xtensor-blas/include/xtensor-blas/xblas.hpp"
#include "../extern/xtensor/include/xtensor/xarray.hpp"
#include "../extern/xtensor-blas/include/xtensor-blas/xlinalg.hpp"
#include "../extern/xtensor-blas/include/xtensor-blas/xlapack.hpp"

// Tests that the Foo::Bar() method does Abc.
TEST(ChainTest, test1) {
  using namespace openbps;
  const std::string input_filepath = "./Xmls/chain_simple.xml";
  Chain chain(read_chain_xml(input_filepath));
  auto react_map = chain.form_reaction();
  EXPECT_EQ(react_map["(n,gamma)"].size(), 4);
}

// Tests that the Foo::Bar() method does Abc.
TEST(ChainTest, test2) {
  using namespace openbps;
  const std::string input_filepath = "./Xmls/chain_simple.xml";
  Chain chain(read_chain_xml(input_filepath));
  auto mapc = chain.form_yield_map();
  EXPECT_EQ(mapc[0.0253][0][0], 6);
  ASSERT_EQ(mapc[0.0253][0][1], 3);
}

// Tests that the Foo::Bar() method does Abc.
TEST(CompositionTest, test3) {
  using namespace openbps;
  pugi::xml_document doc;
  const std::string input_filepath = "./Xmls/reactions.xml";
  auto result = doc.load_file(input_filepath.c_str());
  int res {0};
  pugi::xml_node root_node = doc.child("compositions");
  for (pugi::xml_node tool : root_node.children("composit")) {
      res++;
  }
  EXPECT_EQ(res, 3);
  
}

// Tests that the Foo::Bar() method does Abc.
TEST(ConfigureTest, test4) {
  using namespace openbps;
  pugi::xml_document doc;
  const std::string input_filepath = "./configure.xml";
  auto result = doc.load_file(input_filepath.c_str());
  // Get root element
  pugi::xml_node root = doc.document_element();
  int some = std::stoi(get_node_value(root, "numbers"));
  EXPECT_EQ(some, 1);
  
}

TEST(xlapack, single_element_inverse)
    {
        using namespace xt;
        xtensor<double, 2> a = xt::ones<double>({1, 1});
        auto res = linalg::inv(a);
        EXPECT_EQ(res(), 1.);
    }

TEST(xlapack, test5)
    {
        using namespace xt;
        xarray<double> a = {{ 2, 1, 1},
                            {-1, 1,-1},
                            { 1, 2, 3}};

        xarray<double> vec = {2, 3, -10};
        xarray<double> expected = {3, 1, -5};

        auto res = linalg::solve(a, vec);
        EXPECT_EQ(expected, res);

        vec.reshape({3, 1});
        expected.reshape({3, 1});
        auto res2 = linalg::solve(a, vec);
        EXPECT_EQ(expected, res2);

        xarray<double> vec2 = {6, 2, -10};
        vec2.reshape({3, 1});

        auto res3 = linalg::solve(a, concatenate(xtuple(vec, vec2 * 3), 1));
        xarray<double> expected3 = {{ 3, 16},
                                    { 1,  4},
                                    {-5,-18}};
        EXPECT_EQ(expected3, res3);
    }

