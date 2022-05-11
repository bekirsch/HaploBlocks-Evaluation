#include "utils.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <chrono>

// extracts physical SNP positions
void extract_positions(const char *vcf_path, const char *out_path)
{
	std::ifstream vcf_file;
	std::ofstream out_file;
	
	vcf_file.open(vcf_path);
	std::string line;
	std::vector<std::string> positions;
	while (std::getline(vcf_file, line))
	{
		if (line.front() == '#')
			continue;

		auto split_line = splitString(line, '\t', 2600);
		// skip everything that is not a SNP
		if (split_line[3].size() != 1 || split_line[4].size() != 1)
			continue;
		
		positions.push_back(split_line[1]);
	}
	vcf_file.close();
	
	// clear output file
	out_file.open(out_path, std::ofstream::out | std::ofstream::trunc);
	out_file.close();
	out_file.open(out_path, std::ofstream::out | std::ofstream::app);
	for (auto pos : positions)
		out_file << pos << "\n";
	out_file.close();
	std::cout << positions.size() << " positions written." << std::endl;
}

int main(int argc, char *argv[])
{
	// parse arguments
	if (has_cmd_option(argv, argv+argc, "-h") ||
		!has_cmd_option(argv, argv+argc, "-i") ||
		!has_cmd_option(argv, argv+argc, "-o"))
	{
		std::cout << "Call with <prog> -i vcf_file -o out_file [-h]" << std::endl;
		return EXIT_FAILURE;
	}
	auto in_path = get_cmd_option(argv, argv+argc, "-i");
	auto out_path = get_cmd_option(argv, argv+argc, "-o");

	auto begin = std::chrono::steady_clock::now();
	extract_positions(in_path, out_path);
	auto end = std::chrono::steady_clock::now();
	std::cout << "Parsing took (s) " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << std::endl;
	return EXIT_SUCCESS;
}
