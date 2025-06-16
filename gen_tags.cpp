// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#include <omp.h>
#include <cstring>
#include <ctime>
#include <timer.h>
#include "log.h"

#include "utils.h"

#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

void get_point_num_dim(char *base_data_file, size_t &num_points, size_t &dim) {
  std::ifstream reader(base_data_file, std::ios::binary | std::ios::ate);
  reader.seekg(0);
  int npts_i32, dim_i32;
  reader.read((char *) &npts_i32, sizeof(int));
  reader.read((char *) &dim_i32, sizeof(int));
  num_points = npts_i32;
  dim = dim_i32;
}

using TagT = uint32_t;
template<class T>
void gen_tags(char *base_data_file, char *index, bool sector_aligned) {
  size_t num_points, dim;
  get_point_num_dim(base_data_file, num_points, dim);
  LOG(INFO) << "Loaded " << num_points << " points of dim " << dim;
  std::string tag_file = std::string(index) + "_disk.index.tags";
  TagT *tags = new TagT[num_points];
  for (size_t i = 0; i < num_points; i++) {
    tags[i] = i;
  }
  LOG(INFO) << "Saving tags to " << tag_file;
  if (sector_aligned) {
    pipeann::save_bin_sector_aligned<TagT>((tag_file + ".aligned").c_str(), tags, num_points, 1);
  } else {
    pipeann::save_bin<TagT>(tag_file.c_str(), tags, num_points, 1);
  }
}

int main(int argc, char **argv) {
  if (argc < 5) {
    std::cout << "Correct usage: " << argv[0]
              << " <type[int8/uint8/float]> <base_data_file> <index_file_prefix> <sector_aligned>" << std::endl;
    exit(-1);
  }

  int arg_no = 1;
  char *type = argv[arg_no++];
  char *base_data_file = argv[arg_no++];
  char *index_file_prefix = argv[arg_no++];
  bool sector_aligned = (bool) atoi(argv[arg_no++]);

  if (strcmp(type, "int8") == 0) {
    gen_tags<int8_t>(base_data_file, index_file_prefix, sector_aligned);
  } else if (strcmp(type, "uint8") == 0) {
    gen_tags<uint8_t>(base_data_file, index_file_prefix, sector_aligned);
  } else if (strcmp(type, "float") == 0) {
    gen_tags<float>(base_data_file, index_file_prefix, sector_aligned);
  } else {
    std::cout << "Unknown type: " << type << std::endl;
    exit(-1);
  }
}