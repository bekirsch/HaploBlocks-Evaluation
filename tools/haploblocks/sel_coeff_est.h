#pragma once
#include <string>
#include <tuple>

std::tuple<unsigned int, unsigned int> estimate_selection_coeff(
	const std::string in_path, const std::string out_path,
	const std::string recomb_map_path, const std::string pos_list_path,
	const std::string lookup_path, unsigned int chr_number = 0);
