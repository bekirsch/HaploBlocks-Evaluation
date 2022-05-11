#include "../utils.h"
#include "../pbwt.h"
#include <string>
#include <iostream>

class Config
{
public:
	Config(){};
	
	// output folder
	std::string out_file = "";
	
	// location of the vcf file
	std::string block_path = "";
};

bool parse_command_args(Config &config, int argc, char* argv[])
{
	if (has_cmd_option(argv, argv+argc, "-h")
		|| !has_cmd_option(argv, argv+argc, "--out_file")
		|| !has_cmd_option(argv, argv+argc, "--block_path"))
	{
		std::cerr << "Mandatory args: --out_file --block_path" << std::endl;
		std::cerr << "Optional arg: --min_weight" << std::endl;
		return false;
	}
	
	config.block_path = get_cmd_option(argv, argv+argc, "--block_path");
	config.out_file = get_cmd_option(argv, argv+argc, "--out_file");
	return true;
}

int main(int argc, char* argv[])
{
	Config config;
	if (!parse_command_args(config, argc, argv))
		return EXIT_FAILURE;
	
	calc_haploblocks(config.block_path, config.out_file);
	return EXIT_SUCCESS;
}
