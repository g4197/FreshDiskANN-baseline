// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#include "v2/merge_insert.h"
#include <chrono>
#include <cstring>
#include <iomanip>
#include <omp.h>
#include <pq_flash_index.h>
#include <string.h>
#include <time.h>

#include "log.h"
#include "aux_utils.h"
#include "timer.h"
#include "utils.h"

#ifndef _WINDOWS
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include "linux_aligned_file_reader.h"
#else
#ifdef USE_BING_INFRA
#include "bing_aligned_file_reader.h"
#else
#include "windows_aligned_file_reader.h"
#endif
#endif

#define WARMUP false

void print_stats(std::string category, std::vector<float> percentiles,
                 std::vector<float> results) {
  diskann::cout << std::setw(20) << category << ": " << std::flush;
  for (uint32_t s = 0; s < percentiles.size(); s++) {
    diskann::cout << std::setw(8) << percentiles[s] << "%";
  }
  diskann::cout << std::endl;
  diskann::cout << std::setw(22) << " " << std::flush;
  for (uint32_t s = 0; s < percentiles.size(); s++) {
    diskann::cout << std::setw(9) << results[s];
  }
  diskann::cout << std::endl;
}

template<typename T>
int search_disk_index(int argc, char** argv, diskann::Parameters& paras,
                      diskann::Distance<T>* dist_cmp) {
  // load query bin
  T*                query = nullptr;
  unsigned*         gt_ids = nullptr;
  float*            gt_dists = nullptr;
  uint32_t*         tags = nullptr;
  size_t            query_num, query_dim, query_aligned_dim, gt_num, gt_dim;
  std::vector<_u64> Lvec;

  // bool tags_flag = true;

  int         index = 2;
  std::string index_prefix_path(argv[index++]);
  std::string warmup_query_file = index_prefix_path + "_sample_data.bin";
  int         merge_interval = std::atoi(argv[index++]);
  bool        single_file_index = false;
  _u64        num_nodes_to_cache = std::atoi(argv[index++]);
  _u32        num_threads = std::atoi(argv[index++]);
  _u32        beamwidth = std::atoi(argv[index++]);
  std::string data_bin(argv[index++]);
  std::string query_bin(argv[index++]);
  std::string truthset_bin_prefix(argv[index++]);
  _u64        recall_at = std::atoi(argv[index++]);
  std::string result_output_prefix(argv[index++]);
  std::string dist_metric(argv[index++]);
  bool        use_page_search = std::atoi(argv[index++]) != 0;
  _u32        mem_L = std::atoi(argv[index++]);
  float       ratio = std::atof(argv[index++]);
  int         rounds = std::atoi(argv[index++]);

  diskann::Metric m =
      dist_metric == "cosine" ? diskann::Metric::COSINE : diskann::Metric::L2;
  if (dist_metric != "l2" && m == diskann::Metric::L2) {
    diskann::cout << "Unknown distance metric: " << dist_metric
                  << ". Using default(L2) instead." << std::endl;
  }

  std::string disk_index_tag_file = index_prefix_path + "_disk.index.tags";

  bool calc_recall_flag = false;

  for (int ctr = index; ctr < argc; ctr++) {
    _u64 curL = std::atoi(argv[ctr]);
    if (curL >= recall_at)
      Lvec.push_back(curL);
  }

  if (Lvec.size() == 0) {
    diskann::cout
        << "No valid Lsearch found. Lsearch must be at least recall_at"
        << std::endl;
    return -1;
  }

  diskann::cout << "Search parameters: #threads: " << num_threads << ", ";
  if (beamwidth <= 0)
    diskann::cout << "beamwidth to be optimized for each L value" << std::endl;
  else
    diskann::cout << " beamwidth: " << beamwidth << std::endl;

  diskann::load_aligned_bin<T>(query_bin, query, query_num, query_dim,
                               query_aligned_dim);

  std::shared_ptr<AlignedFileReader> reader = nullptr;
#ifdef _WINDOWS
#ifndef USE_BING_INFRA
  reader.reset(new WindowsAlignedFileReader());
#else
  reader.reset(new diskann::BingAlignedFileReader());
#endif
#else
  reader.reset(new LinuxAlignedFileReader());
#endif

  paras.Set<_u32>("nodes_to_cache", num_nodes_to_cache);
  paras.Set<_u32>("num_search_threads", num_threads);
  paras.Set<_u32>("beamwidth", beamwidth);
  diskann::MergeInsert<T> sync_index(
      paras, query_dim, index_prefix_path + "_mem", index_prefix_path,
      index_prefix_path + "_merge", dist_cmp, m, single_file_index,
      index_prefix_path);

  T*     data;
  size_t data_npts, data_dim;
  diskann::load_bin<T>(data_bin, data, data_npts, data_dim);

  std::string recall_string = "Recall@" + std::to_string(recall_at);
  diskann::cout << std::setw(6) << "L" << std::setw(12) << "Beamwidth"
                << std::setw(16) << "QPS" << std::setw(16) << "Mean Latency"
                << std::setw(16) << "99.9 Latency" << std::setw(16)
                << "Mean IOs" << std::setw(16) << "CPU (s)" << std::setw(16)
                << "IO (us)";
  if (calc_recall_flag) {
    diskann::cout << std::setw(16) << recall_string << std::endl;
  } else {
    diskann::cout << std::endl;
  }
  diskann::cout
      << "==============================================================="
         "==========================================="
      << std::endl;

  calc_recall_flag = true;
  uint64_t index_npts = sync_index._disk_index->num_points;
  uint64_t vector_per_round = index_npts * ratio;

  auto st = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < 100; ++i) {
    uint64_t r_start = i * vector_per_round;
    uint64_t r_end = r_start + index_npts;

#pragma omp parallel for num_threads(num_threads)
    for (int j = r_end; j < r_end + vector_per_round; ++j) {
      uint32_t vtag = j;
      T*       v = data + vtag * data_dim;
      sync_index.insert(v, vtag);
      // sync_index.insert_in_place(v, vtag);
    }
    if (i % merge_interval == merge_interval - 1) {
      sync_index.final_merge();
    }
    auto ed = std::chrono::high_resolution_clock::now();
    LOG(INFO) << "Inserted " << (i + 1) * vector_per_round << ", throughput: "
              << 1000 * (double) (i + 1) * vector_per_round /
                     std::chrono::duration_cast<std::chrono::milliseconds>(ed -
                                                                           st)
                         .count()
              << " vec/s";
  }
  diskann::aligned_free(query);
  delete[] gt_ids;
  delete[] gt_dists;
  return 0;
}

int main(int argc, char** argv) {
  if (argc < 15) {
    // tags == 1!
    diskann::cout
        << "Usage: " << argv[0]
        << " <index_type (float/int8/uint8)>  <index_prefix_path>"
           " <merge_interval>"
           " <num_nodes_to_cache>  <num_threads>  <beamwidth (use 0 to "
           "optimize internally)> "
           " <data.bin> <query_file.bin>  <truthset prefix> "
           " <K>  <result_output_prefix> <similarity (cosine/l2)> "
           " <use_page_search(0/1)> <mem_L> <ratio> <rounds> <L1> [L2] etc.  "
           "See README for "
           "more information on parameters."
        << std::endl;
    exit(-1);
  }

  diskann::Parameters paras;
  paras.Set<unsigned>("L_mem", 128);
  paras.Set<unsigned>("R_mem", 48);
  paras.Set<float>("alpha_mem", 1.2);
  paras.Set<unsigned>("L_disk", 128);
  paras.Set<unsigned>("R_disk", 48);
  paras.Set<float>("alpha_disk", 1.2);
  paras.Set<unsigned>("C", 75);

  if (std::string(argv[1]) == std::string("float")) {
    diskann::DistanceL2 dist_cmp;
    search_disk_index<float>(argc, argv, paras, &dist_cmp);
  } else if (std::string(argv[1]) == std::string("int8")) {
    diskann::DistanceL2Int8 dist_cmp;
    search_disk_index<int8_t>(argc, argv, paras, &dist_cmp);
  } else if (std::string(argv[1]) == std::string("uint8")) {
    diskann::DistanceL2UInt8 dist_cmp;
    search_disk_index<uint8_t>(argc, argv, paras, &dist_cmp);
  } else {
    diskann::cout << "Unsupported index type. Use float or int8 or uint8"
                  << std::endl;
  }
}
