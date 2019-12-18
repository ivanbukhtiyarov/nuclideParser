#include "nuclide_class.h"
#include "configure.h"
#include "gtest/gtest.h"
#include "reactions.h"
#include "../extern/pugiData/pugixml.h"


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
  EXPECT_EQ(some, 30);
  
}

