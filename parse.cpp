#include "nuclide_class.h"

std::vector<std::string> split(const std::string& s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<double> splitAtof(const std::string& s, char delimiter)
{
    std::vector<double> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {

        tokens.push_back(atof(token.c_str()));
    }
    return tokens;
}

t_decay parse_decay(pugi::xml_node node)
{
    t_decay temp;

    temp.type = node.attribute("type").value();
    temp.target = node.attribute("target").value();
    temp.branching_ratio = atof(node.attribute("branching_ratio").value());
    return temp;
}

t_reaction parse_reaction(pugi::xml_node node)
{
    t_reaction temp;

    temp.type = node.attribute("type").value();
    temp.q = atof(node.attribute("Q").value());
    temp.target = node.attribute("target").value();
    return temp;
}

t_yield parse_yield(pugi::xml_node node)
{
    t_yield                  temp;
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

t_nfy parse_nfy(pugi::xml_node node)
{
    t_nfy           temp;

    temp.energies = splitAtof(node.child("energies").child_value(), ' ');
    for (pugi::xml_node tool = node.child("fission_yields"); tool; tool = tool.next_sibling("fission_yields"))
    {
        temp.yield_arr.push_back(parse_yield(tool));
    }
    return temp;
}

t_nuclide parse_nuclide(pugi::xml_node node)
{
    t_nuclide       temp;

    temp.name = node.attribute("name").value();
    temp.decay_modes = atoi(node.attribute("decay_modes").value());
    temp.reactions = atoi(node.attribute("reactions").value());
    temp.half_life = atof(node.attribute("half_life").value());
    if (temp.decay_modes > 0)
    {
        temp.decay_energy = atof(node.attribute("decay_energy").value());
        for (pugi::xml_node tool = node.child("decay"); tool; tool = tool.next_sibling("decay"))
        {
            temp.decay_arr.push_back(parse_decay(tool));
        }
    }
    else
    {
        temp.decay_energy = 0;
    }

    if (temp.reactions > 0)
    {
        for (pugi::xml_node tool = node.child("reaction"); tool; tool = tool.next_sibling("reaction"))
        {
            temp.reaction_arr.push_back(parse_reaction(tool));
        }
    }

    if(node.child("neutron_fission_yields"))
    {
        temp.nfy = parse_nfy(node.child("neutron_fission_yields"));
    }
    return temp;
}

t_chain parse_chain(pugi::xml_node node)
{
    t_chain temp;
    t_nuclide temp_nuclide;

    for (pugi::xml_node tool = node.child("nuclide"); tool; tool = tool.next_sibling("nuclide"))
    {
        temp_nuclide = parse_nuclide(tool);
        temp.nuclides.insert({temp_nuclide.name, temp_nuclide});
    }
    return temp;
}