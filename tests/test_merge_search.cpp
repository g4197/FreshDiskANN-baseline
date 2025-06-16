// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

#include "v2/index_merger.h"
#include "v2/merge_insert.h"

#include <atomic>
#include <cstring>
#include <iomanip>
#include <omp.h>
#include <pq_flash_index.h>
#include <set>
#include <string.h>
#include <time.h>

#include "log.h"
#include "aux_utils.h"
#include "index.h"
#include "math_utils.h"
#include "memory_mapper.h"
#include "partition_and_pq.h"
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
  diskann::cerr << std::setw(20) << category << ": " << std::flush;
  for (uint32_t s = 0; s < percentiles.size(); s++) {
    diskann::cerr << std::setw(8) << percentiles[s] << "%";
  }
  diskann::cerr << std::endl;
  diskann::cerr << std::setw(22) << " " << std::flush;
  for (uint32_t s = 0; s < percentiles.size(); s++) {
    diskann::cerr << std::setw(9) << results[s];
  }
  diskann::cerr << std::endl;
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
  bool        single_file_index = std::atoi(argv[index++]) != 0;
  _u64        num_nodes_to_cache = std::atoi(argv[index++]);
  _u32        num_threads = std::atoi(argv[index++]);
  _u32        beamwidth = std::atoi(argv[index++]);
  std::string data_bin(argv[index++]);
  std::string query_bin(argv[index++]);
  std::string truthset_bin(argv[index++]);
  _u64        recall_at = std::atoi(argv[index++]);
  std::string result_output_prefix(argv[index++]);
  std::string dist_metric(argv[index++]);
  bool        use_page_search = std::atoi(argv[index++]) != 0;
  _u32        mem_L = std::atoi(argv[index++]);
  int         vector_per_round = std::atoi(argv[index++]);
  int         rounds = std::atoi(argv[index++]);

  diskann::Metric m =
      dist_metric == "cosine" ? diskann::Metric::COSINE : diskann::Metric::L2;
  if (dist_metric != "l2" && m == diskann::Metric::L2) {
    diskann::cerr << "Unknown distance metric: " << dist_metric
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
    diskann::cerr
        << "No valid Lsearch found. Lsearch must be at least recall_at"
        << std::endl;
    return -1;
  }

  diskann::cerr << "Search parameters: #threads: " << num_threads << ", ";
  if (beamwidth <= 0)
    diskann::cerr << "beamwidth to be optimized for each L value" << std::endl;
  else
    diskann::cerr << " beamwidth: " << beamwidth << std::endl;

  diskann::load_aligned_bin<T>(query_bin, query, query_num, query_dim,
                               query_aligned_dim);

  if (file_exists(truthset_bin)) {
    diskann::load_truthset(truthset_bin, gt_ids, gt_dists, gt_num, gt_dim,
                           &tags);
    if (gt_num != query_num) {
      diskann::cerr
          << "Error. Mismatch in number of queries and ground truth data"
          << std::endl;
    }
    calc_recall_flag = true;
    LOG(INFO) << "Loaded ground truth num " << gt_num << " with dim " << gt_dim;
  }

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
  diskann::cerr << std::setw(6) << "L" << std::setw(12) << "Beamwidth"
                << std::setw(16) << "QPS" << std::setw(16) << "Mean Latency"
                << std::setw(16) << "99.9 Latency" << std::setw(16)
                << "Mean IOs" << std::setw(16) << "CPU (s)" << std::setw(16)
                << "IO (us)";
  if (calc_recall_flag) {
    diskann::cerr << std::setw(16) << recall_string << std::endl;
  } else {
    diskann::cerr << std::endl;
  }
  diskann::cerr
      << "==============================================================="
         "==========================================="
      << std::endl;

  for (int i = 0; i < rounds; ++i) {
    // search.
    std::vector<std::vector<uint32_t>> query_result_ids(Lvec.size());
    std::vector<std::vector<uint32_t>> query_result_tags(Lvec.size());
    std::vector<std::vector<float>>    query_result_dists(Lvec.size());
    for (uint32_t test_id = 0; test_id < Lvec.size(); test_id++) {
      _u64     L = Lvec[test_id];
      uint32_t optimized_beamwidth = beamwidth;

      query_result_ids[test_id].resize(recall_at * query_num);
      query_result_dists[test_id].resize(recall_at * query_num);
      query_result_tags[test_id].resize(recall_at * query_num);

      diskann::QueryStats* stats = new diskann::QueryStats[query_num];

      std::vector<uint64_t> query_result_tags_64(recall_at * query_num);
      std::vector<uint32_t> query_result_tags_32(recall_at * query_num);
      auto                  s = std::chrono::high_resolution_clock::now();

#pragma omp parallel for schedule(dynamic, 1)
      for (_s64 i = 0; i < (int64_t) query_num; i++) {
        sync_index.search_sync(
            query + (i * query_aligned_dim), (uint64_t) recall_at, (uint64_t) L,
            query_result_tags_32.data() + (i * recall_at),
            query_result_dists[test_id].data() + (i * recall_at), stats + i);
      }

      auto e = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> diff = e - s;
      float                         qps =
          (float) ((1.0 * (double) query_num) / (1.0 * (double) diff.count()));

      diskann::convert_types<uint32_t, uint32_t>(
          query_result_tags_32.data(), query_result_tags[test_id].data(),
          (size_t) query_num, (size_t) recall_at);

      float mean_latency = (float) diskann::get_mean_stats(
          stats, query_num,
          [](const diskann::QueryStats& stats) { return stats.total_us; });

      float latency_999 = (float) diskann::get_percentile_stats(
          stats, query_num, 0.999f,
          [](const diskann::QueryStats& stats) { return stats.total_us; });

      float mean_ios = (float) diskann::get_mean_stats(
          stats, query_num,
          [](const diskann::QueryStats& stats) { return stats.n_ios; });

      float mean_cpuus = (float) diskann::get_mean_stats(
          stats, query_num,
          [](const diskann::QueryStats& stats) { return stats.cpu_us; });

      float mean_ious = (float) diskann::get_mean_stats(
          stats, query_num,
          [](const diskann::QueryStats& stats) { return stats.io_us; });
      delete[] stats;

      float recall = 0;
      if (calc_recall_flag) {
        recall = (float) diskann::calculate_recall(
            (_u32) query_num, gt_ids, gt_dists, (_u32) gt_dim,
            query_result_tags[test_id].data(), (_u32) recall_at,
            (_u32) recall_at);
      }

      diskann::cerr << std::setw(6) << L << std::setw(12) << optimized_beamwidth
                    << std::setw(16) << qps << std::setw(16) << mean_latency
                    << std::setw(16) << latency_999 << std::setw(16) << mean_ios
                    << std::setw(16) << mean_cpuus << std::setw(16)
                    << mean_ious;
      if (calc_recall_flag) {
        diskann::cerr << std::setw(16) << recall << std::endl;
      }
      // randomly select several IDs for removal and insertion.
      // LOG(INFO) << "Remove and insert round " << i;
      // for (int j = i * vector_per_round; j < (i + 1) * vector_per_round; ++j)
      // {
      //   sync_index.lazy_delete(vectors_to_update[j].first);
      // }

      auto st = std::chrono::high_resolution_clock::now();
      for (int j = i * vector_per_round; j < (i + 1) * vector_per_round; ++j) {
        uint32_t vtag = j % data_npts;
        sync_index.lazy_delete(vtag);
        // sync_index.delete_in_place(vtag);
      }
      for (int j = i * vector_per_round; j < (i + 1) * vector_per_round; ++j) {
        uint32_t vtag = j % data_npts;
        T*       v = data + vtag * data_dim;
        sync_index.insert(v, vtag);
        // sync_index.insert_in_place(v, vtag);
      }
      sync_index.final_merge();
      // T dummy[query_dim];
      // sync_index.insert(dummy, npts + 1);
      // sync_index.final_merge();
      // sync_index.final_merge();
    }
  }
  diskann::aligned_free(query);
  delete[] gt_ids;
  delete[] gt_dists;
  return 0;
}

int main(int argc, char** argv) {
  if (argc < 15) {
    // tags == 1!
    diskann::cerr
        << "Usage: " << argv[0]
        << " <index_type (float/int8/uint8)>  <index_prefix_path>"
           " <single_file_index(0/1)>"
           " <num_nodes_to_cache>  <num_threads>  <beamwidth (use 0 to "
           "optimize internally)> "
           " <data.bin> <query_file.bin>  <truthset.bin (use \"null\" for "
           "none)> "
           " <K>  <result_output_prefix> <similarity (cosine/l2)> "
           " <use_page_search(0/1)> <mem_L> <vector_per_round> <rounds> <L1> "
           "[L2] etc.  See README for "
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
    diskann::cerr << "Unsupported index type. Use float or int8 or uint8"
                  << std::endl;
  }
}
