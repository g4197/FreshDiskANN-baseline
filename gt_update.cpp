// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#include <omp.h>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <timer.h>
#include "log.h"
#include "utils.h"

#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

constexpr double ratio = 0.01;
void check_and_gen(char *gt_path, uint64_t tot_npts, uint64_t batch_npts, uint64_t target_topk, char *target_dir,
                   bool insert_only) {
  // st: [0, tot_npts / 2]; ed: [tot_npts / 2, tot_npts]
  // expect: gt is more than top 100.
  uint32_t *data = nullptr;
  size_t nq, dim;
  pipeann::load_bin<uint32_t>(gt_path, data, nq, dim);
  auto data_idx = [&](size_t x, size_t y) { return data[x * dim + y]; };
  LOG(INFO) << "Loaded " << nq << " points with dim " << dim << " from " << gt_path;
  LOG(INFO) << "Checking if gt is more than top " << target_topk << " for each query.";
  uint64_t step = ratio * batch_npts;
  int success_cnt = 0;
  for (uint64_t st = 0; st < (tot_npts - batch_npts); st += step) {
    bool success = true;
    auto real_st = insert_only ? 0 : st;
    uint64_t ed = st + batch_npts;
    LOG(INFO) << "Checking range [" << real_st << ", " << ed << ")";
    std::vector<uint32_t> cur_gt;
    for (uint64_t i = 0; i < nq; ++i) {
      int cnt = 0;
      for (uint64_t j = 0; j < dim; ++j) {
        if (data_idx(i, j) >= real_st && data_idx(i, j) < ed) {
          cur_gt.push_back(data_idx(i, j));
          ++cnt;
        }
        if ((uint64_t) cnt >= target_topk) {
          break;
        }
      }
      if ((uint64_t) cnt < target_topk) {
        LOG(FATAL) << "Query " << i << " has less than " << target_topk << " gt points in range [" << st << ", " << ed
                   << ")"
                   << " " << cnt;
        success = false;
      }
    }
    success_cnt += success;
    if (success) {
      std::string target_path = std::string(target_dir) + "/gt_" + std::to_string(st) + ".bin";
      LOG(INFO) << "Writing to " << target_path;
      pipeann::save_bin<uint32_t>(target_path.c_str(), cur_gt.data(), nq, target_topk);
    }
  }
  LOG(INFO) << "Success rate: " << success_cnt << " / " << ((tot_npts - batch_npts) / step);
}

int main(int argc, char **argv) {
  if (argc < 6) {
    std::cout << "Correct usage: " << argv[0]
              << " <file> <tot_npts> <batch_npts> <target_topk> <target_dir> <insert_only>" << std::endl;
    exit(-1);
  }

  int arg_no = 1;
  char *gt_path = argv[arg_no++];
  uint64_t npts = std::stoull(argv[arg_no++]);
  uint64_t batch_npts = std::stoull(argv[arg_no++]);
  uint64_t target_topk = std::stoull(argv[arg_no++]);
  char *target_dir = argv[arg_no++];
  bool insert_only = std::stoi(argv[arg_no++]);

  // mkdir target_dir
  std::string cmd = "mkdir -p " + std::string(target_dir);
  std::ignore = system(cmd.c_str());
  check_and_gen(gt_path, npts, batch_npts, target_topk, target_dir, insert_only);
}