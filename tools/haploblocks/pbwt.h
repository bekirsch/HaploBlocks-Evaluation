#pragma once
#include <string>

typedef long long int int_all;
int_all calc_haploblocks(std::string block_file, std::string out_path,
    int_all minimal_block_size = 2, int_all buffer_size = 8192);
