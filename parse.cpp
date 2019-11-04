#include "nuclide_class.h"

namespace openbps {

std::string
get_node_value(pugi::xml_node node, const char* name)
{
  // Search for either an attribute or child tag and get the data as a char*.
  const pugi::char_t* value_char;
  if (node.attribute(name)) {
    value_char = node.attribute(name).value();
  } else if (node.child(name)) {
    value_char = node.child_value(name);
  } else {
    
    std::cerr << "Node \"" << name << "\" is not a member of the \""
              << node.name() << "\" XML node";
    
  }
  std::string value {value_char};
  return value;
}

bool
get_node_value_bool(pugi::xml_node node, const char* name)
{
  if (node.attribute(name)) {
    return node.attribute(name).as_bool();
  } else if (node.child(name)) {
    return node.child(name).text().as_bool();
  } else {
    
    std::cerr << "Node \"" << name << "\" is not a member of the \""
              << node.name() << "\" XML node";
    
  }
  return false;
}

std::vector<std::string> split(const std::string& s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<double> splitAtof(const std::string& s, char delimiter)
{
    std::vector<double> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(atof(token.c_str()));
    }
    return tokens;
}


//==============================================================================
// Chain implementation
//==============================================================================


s_decay Chain::parse_decay_(pugi::xml_node node)
{
    s_decay temp;

    temp.type = node.attribute("type").value();
    temp.target = node.attribute("target").value();
    temp.branching_ratio = atof(node.attribute("branching_ratio").value());
    return temp;
}

s_reaction Chain::parse_reaction_(pugi::xml_node node)
{
    s_reaction temp;

    temp.type = node.attribute("type").value();
    temp.q = atof(node.attribute("Q").value());
    temp.target = node.attribute("target").value();
    return temp;
}

s_yield Chain::parse_yield_(pugi::xml_node node)
{
    s_yield                  temp;
    std::vector<std::string> nuclides;
    std::vector<double>      numbers;

    temp.energy = atof(node.attribute("energy").value());
    nuclides = split(node.child("products").child_value(),' ');
    numbers = splitAtof(node.child("data").child_value(),' ');
    for(std::size_t i = 0; i < nuclides.size(); i++)
    {
        temp.product_data.insert({nuclides[i],numbers[i]});
    }
    return temp;
}

s_nfy Chain::parse_nfy_(pugi::xml_node node)
{
    s_nfy temp;

    temp.energies = splitAtof(node.child("energies").child_value(), ' ');
    for (pugi::xml_node tool : node.children("fission_yields")) {
        temp.yield_arr.push_back(parse_yield_(tool));
    }
    return temp;
}

s_nuclide Chain::parse_nuclide_(pugi::xml_node node)
{
    s_nuclide temp;

    temp.name = node.attribute("name").value();
    temp.decay_modes = atoi(node.attribute("decay_modes").value());
    temp.reactions = atoi(node.attribute("reactions").value());
    temp.half_life = atof(node.attribute("half_life").value());
    if (temp.decay_modes > 0) {
        temp.decay_energy = atof(node.attribute("decay_energy").value());
        for (pugi::xml_node tool : node.children("decay")) {
            temp.decay_arr.push_back(parse_decay_(tool));
        }
    }
    else {
        temp.decay_energy = 0;
    }

    if (temp.reactions > 0) {
        for (pugi::xml_node tool : node.children("reaction")) {
            temp.reaction_arr.push_back(parse_reaction_(tool));
        }
    }

    if(node.child("neutron_fission_yields")) {
        temp.nfy = parse_nfy_(node.child("neutron_fission_yields"));
    }
    return temp;
}

Chain::Chain(pugi::xml_node node)
{
    
    for (pugi::xml_node tool : node.children("nuclide")) {
        nuclides.insert({get_node_value(tool, "name"), parse_nuclide_(tool)});
    }
    
}

} //namespace openbps

