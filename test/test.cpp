#include "nuclide_class.h"

#include "gtest/gtest.h"



// Tests that the Foo::Bar() method does Abc.
TEST(ChainTest, test1) {
  using namespace openbps;
  const std::string input_filepath = "./Xmls/chain_simple.xml";
  Chain chain(read_xml(input_filepath));
  auto react_map = chain.form_reaction();
  EXPECT_EQ(react_map["(n,gamma)"].size(), 4);
}

// Tests that the Foo::Bar() method does Abc.
TEST(ChainTest, test2) {
  using namespace openbps;
  const std::string input_filepath = "./Xmls/chain_simple.xml";
  Chain chain(read_xml(input_filepath));
  auto mapc = chain.form_yield_map();
  EXPECT_EQ(mapc["0.0253"][0][0], 6);
  ASSERT_EQ(mapc["0.0253"][0][1], 3);
}

