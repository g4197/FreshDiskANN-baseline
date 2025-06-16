#include "v2/merge_insert.h"

#include "log.h"
#include <index.h>
#include <cstddef>
#include <future>
#include <Neighbor_Tag.h>
#include <numeric>
#include <omp.h>
#include <string.h>
#include <time.h>
#include <timer.h>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <dirent.h>
#include <sys/stat.h>

#include "aux_utils.h"
#include "index.h"
#include "math_utils.h"
#include "partition_and_pq.h"
#include "utils.h"

#ifndef _WINDOWS
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

uint32_t Merge_Size = 0;
#define NUM_INSERT_THREADS 24
#define NUM_DELETE_THREADS 1
#define NUM_SEARCH_THREADS 32
#define InsertQPS 1200
#define DeleteQPS 1200

int            begin_time = 0;
diskann::Timer globalTimer;

// acutually also shows disk size
void ShowMemoryStatus() {
  int current_time = globalTimer.elapsed() / 1.0e6f - begin_time;

  int           tSize = 0, resident = 0, share = 0;
  std::ifstream buffer("/proc/self/statm");
  buffer >> tSize >> resident >> share;
  buffer.close();
  long page_size_kb = sysconf(_SC_PAGE_SIZE) /
                      1024;  // in case x86-64 is configured to use 2MB pages
  double rss = resident * page_size_kb;

  std::cout << "memory current time: " << current_time << " RSS : " << rss
            << " KB" << std::endl;
}

std::string convertFloatToString(const float value, const int precision = 0) {
  std::stringstream stream{};
  stream << std::fixed << std::setprecision(precision) << value;
  return stream.str();
}

std::string GetTruthFileName(std::string& truthFilePrefix, int r_start) {
  std::string fileName(truthFilePrefix);
  fileName = fileName + "/gt_" + std::to_string(r_start) + ".bin";
  LOG(INFO) << "Truth file name: " << fileName;
  return fileName;
}

template<typename T>
inline uint64_t save_bin_test(const std::string& filename, T* id, float* dist,
                              size_t npts, size_t ndims, size_t offset = 0) {
  std::ofstream writer;
  open_file_to_write(writer, filename);

  diskann::cout << "Writing bin: " << filename.c_str() << std::endl;
  writer.seekp(offset, writer.beg);
  int    npts_i32 = (int) npts, ndims_i32 = (int) ndims;
  size_t bytes_written = npts * ndims * sizeof(T) + 2 * sizeof(uint32_t);
  writer.write((char*) &npts_i32, sizeof(int));
  writer.write((char*) &ndims_i32, sizeof(int));
  diskann::cout << "bin: #pts = " << npts << ", #dims = " << ndims
                << ", size = " << bytes_written << "B" << std::endl;

  for (int i = 0; i < npts; i++) {
    for (int j = 0; j < ndims; j++) {
      writer.write((char*) (id + i * ndims + j), sizeof(T));
      writer.write((char*) (dist + i * ndims + j), sizeof(float));
    }
  }
  writer.close();
  diskann::cout << "Finished writing bin." << std::endl;
  return bytes_written;
}

template<typename T, typename TagT>
void sync_search_kernel(T* query, size_t query_num, size_t query_aligned_dim,
                        const int recall_at, _u64 L,
                        diskann::MergeInsert<T, TagT>& sync_index,
                        std::string& truthset_file, bool merged,
                        bool calRecall) {
  unsigned* gt_ids = NULL;
  float*    gt_dists = NULL;
  size_t    gt_num, gt_dim;

  if (!file_exists(truthset_file)) {
    calRecall = false;
  }

  if (calRecall) {
    std::cout << "current truthfile: " << truthset_file << std::endl;
    diskann::load_truthset(truthset_file, gt_ids, gt_dists, gt_num, gt_dim);
  }

  float* query_result_dists = new float[recall_at * query_num];
  TagT*  query_result_tags = new TagT[recall_at * query_num];

  for (_u32 q = 0; q < query_num; q++) {
    for (_u32 r = 0; r < (_u32) recall_at; r++) {
      query_result_tags[q * recall_at + r] = std::numeric_limits<TagT>::max();
      query_result_dists[q * recall_at + r] = std::numeric_limits<float>::max();
    }
  }

  std::vector<double> latency_stats(query_num, 0);
  diskann::QueryStats stats[query_num];
  std::string         recall_string = "Recall@" + std::to_string(recall_at);
  std::cout << std::setw(4) << "Ls" << std::setw(12) << "QPS " << std::setw(18)
            << "Mean Lat" << std::setw(12) << "50 Lat" << std::setw(12)
            << "90 Lat" << std::setw(12) << "95 Lat" << std::setw(12)
            << "99 Lat" << std::setw(12) << "99.9 Lat" << std::setw(12)
            << recall_string << std::setw(12) << "Disk IOs" << std::endl;
  std::cout << "==============================================================="
               "==============="
            << std::endl;
  auto s = std::chrono::high_resolution_clock::now();
#pragma omp parallel for num_threads(NUM_SEARCH_THREADS) schedule(dynamic)
  for (int64_t i = 0; i < (int64_t) query_num; i++) {
    auto qs = std::chrono::high_resolution_clock::now();
    // stats[i].n_current_used = 8;
    sync_index.search_sync(query + i * query_aligned_dim, recall_at, L,
                           query_result_tags + i * recall_at,
                           query_result_dists + i * recall_at, stats + i);

    auto qe = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = qe - qs;
    latency_stats[i] = diff.count() * 1000;
    std::this_thread::sleep_for(std::chrono::milliseconds(2));
  }
  auto e = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> diff = e - s;
  float                         qps = (query_num / diff.count());
  float                         recall = 0;

  int current_time = globalTimer.elapsed() / 1.0e6f - begin_time;
  if (calRecall) {
    // for (int i = 0; i < recall_at; ++i) {
    //   LOG(INFO) << query_result_tags[i] << " " << gt_ids[i];
    // }
    recall = diskann::calculate_recall(query_num, gt_ids, gt_dists, gt_dim,
                                       query_result_tags, recall_at, recall_at);
    delete[] gt_ids;
    gt_ids = nullptr;
  }

  std::cout << "search current time: " << current_time << std::endl;

  float mean_ios = (float) diskann::get_mean_stats(
      stats, query_num,
      [](const diskann::QueryStats& stats) { return stats.n_ios; });

  std::sort(latency_stats.begin(), latency_stats.end());
  std::cout << std::setw(4) << L << std::setw(12) << qps << std::setw(18)
            << ((float) std::accumulate(latency_stats.begin(),
                                        latency_stats.end(), 0)) /
                   (float) query_num
            << std::setw(12)
            << (float) latency_stats[(_u64) (0.50 * ((double) query_num))]
            << std::setw(12)
            << (float) latency_stats[(_u64) (0.90 * ((double) query_num))]
            << std::setw(12)
            << (float) latency_stats[(_u64) (0.95 * ((double) query_num))]
            << std::setw(12)
            << (float) latency_stats[(_u64) (0.99 * ((double) query_num))]
            << std::setw(12)
            << (float) latency_stats[(_u64) (0.999 * ((double) query_num))]
            << std::setw(12) << recall << std::setw(12) << mean_ios
            << std::endl;

  delete[] query_result_dists;
  delete[] query_result_tags;
}

template<typename T, typename TagT>
void merge_kernel(diskann::MergeInsert<T, TagT>& sync_index,
                  std::string&                   save_path) {
  // sync_index.dummy_merge();
  //   ShowMemoryStatus();
  sync_index.final_merge();
}

template<typename T, typename TagT>
void deletion_kernel(T* data_load, diskann::MergeInsert<T, TagT>& sync_index,
                     std::vector<TagT>& delete_vec, size_t aligned_dim) {
  diskann::Timer      timer;
  size_t              npts = delete_vec.size();
  std::vector<double> delete_latencies(npts, 0);
  std::cout << "Begin Delete" << std::endl;
  std::atomic_size_t success(0);
#pragma omp parallel for num_threads(NUM_DELETE_THREADS)
  for (_s64 i = 0; i < (_s64) delete_vec.size(); i++) {
    diskann::Timer      delete_timer;
    diskann::QueryStats stats;
    sync_index.lazy_delete(delete_vec[i]);
    success++;
    delete_latencies[i] = ((double) delete_timer.elapsed());

    if (success == DeleteQPS) {
      std::this_thread::sleep_for(std::chrono::seconds(1));
      success = 0;
    }
  }
  float time_secs = timer.elapsed() / 1.0e6f;
  std::sort(delete_latencies.begin(), delete_latencies.end());
  // std::cout << "Deleted " << delete_vec.size() << " / " << success.load() <<
  // " points in " << time_secs << "s"
  //           << std::endl;
  std::cout << "Mem index deletion time : " << timer.elapsed() / 1000 << " ms"
            << std::endl
            << "10th percentile deletion time : "
            << delete_latencies[(size_t) (0.10 * ((double) npts))]
            << " microsec" << std::endl
            << "50th percentile deletion time : "
            << delete_latencies[(size_t) (0.5 * ((double) npts))] << " microsec"
            << std::endl
            << "90th percentile deletion time : "
            << delete_latencies[(size_t) (0.90 * ((double) npts))]
            << " microsec" << std::endl
            << "99th percentile deletion time : "
            << delete_latencies[(size_t) (0.99 * ((double) npts))]
            << " microsec" << std::endl
            << "99.9th percentile deletion time : "
            << delete_latencies[(size_t) (0.999 * ((double) npts))]
            << " microsec" << std::endl;
}

template<typename T, typename TagT>
void insertion_kernel(T* data_load, diskann::MergeInsert<T, TagT>& sync_index,
                      std::vector<TagT>& insert_vec, size_t aligned_dim) {
  diskann::Timer      timer;
  size_t              npts = insert_vec.size();
  std::vector<double> insert_latencies(npts, 0);
  std::cout << "Begin Insert" << std::endl;
  std::atomic_size_t success(0);

  auto st = std::chrono::high_resolution_clock::now();

  auto getQPS = [&]() {
    auto cur = std::chrono::high_resolution_clock::now();
    auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(cur - st);
    if (dur.count() == 0) {
      return 0ul;
    }
    return (1000 * success) / dur.count();
  };

#pragma omp parallel for num_threads(2)
  for (_s64 i = 0; i < (_s64) insert_vec.size(); i++) {
    diskann::Timer insert_timer;
    sync_index.insert(data_load + aligned_dim * i, insert_vec[i]);
    ++success;
    insert_latencies[i] = ((double) insert_timer.elapsed());
    while (getQPS() >= InsertQPS) {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
  }
  float time_secs = timer.elapsed() / 1.0e6f;
  std::sort(insert_latencies.begin(), insert_latencies.end());
  std::cout << "Inserted " << insert_vec.size() << " points in " << time_secs
            << "s" << std::endl;
  std::cout << "Mem index insertion time : " << timer.elapsed() / 1000 << " ms"
            << std::endl
            << "10th percentile insertion time : "
            << insert_latencies[(size_t) (0.10 * ((double) npts))]
            << " microsec" << std::endl
            << "50th percentile insertion time : "
            << insert_latencies[(size_t) (0.5 * ((double) npts))] << " microsec"
            << std::endl
            << "90th percentile insertion time : "
            << insert_latencies[(size_t) (0.90 * ((double) npts))]
            << " microsec" << std::endl
            << "99th percentile insertion time : "
            << insert_latencies[(size_t) (0.99 * ((double) npts))]
            << " microsec" << std::endl
            << "99.9th percentile insertion time : "
            << insert_latencies[(size_t) (0.999 * ((double) npts))]
            << " microsec" << std::endl;
}

template<typename T, typename TagT = uint32_t>
void get_trace(std::string data_bin, uint64_t l_start, uint64_t r_start,
               uint64_t n, std::vector<TagT>& delete_tags,
               std::vector<TagT>& insert_tags, std::vector<T>& data_load) {
  for (uint64_t i = l_start; i < l_start + n; ++i) {
    delete_tags.push_back(i);
  }

  for (uint64_t i = r_start; i < r_start + n; ++i) {
    insert_tags.push_back(i);
  }

  // load data, load n vecs from r_start.
  int           npts_i32, dim_i32;
  std::ifstream reader(data_bin, std::ios::binary | std::ios::ate);
  reader.seekg(0, reader.beg);
  reader.read((char*) &npts_i32, sizeof(int));
  reader.read((char*) &dim_i32, sizeof(int));

  size_t data_dim = dim_i32;
  data_load.resize(n * data_dim);
  reader.seekg(2 * sizeof(int) + r_start * data_dim * sizeof(T), reader.beg);
  reader.read((char*) data_load.data(), sizeof(T) * n * data_dim);
}

template<typename T, typename TagT>
void update(const std::string& data_bin, const unsigned L_disk,
            const unsigned R_disk, const float alpha_disk, int step,
            const size_t base_num, const unsigned nodes_to_cache,
            std::string& save_path, const std::string& query_file,
            std::string& truthset_file, const int recall_at,
            std::vector<_u64> Lsearch, const unsigned beam_width,
            diskann::Distance<T>* dist_cmp, int ckpt = 0) {
  diskann::Parameters paras;
  paras.Set<unsigned>("L_mem", 128);
  paras.Set<unsigned>("R_mem", 64);
  paras.Set<float>("alpha_mem", 1.2);
  paras.Set<unsigned>("L_disk", L_disk);
  paras.Set<unsigned>("R_disk", R_disk);
  paras.Set<float>("alpha_disk", alpha_disk);
  paras.Set<unsigned>("C", 384);
  paras.Set<unsigned>("beamwidth", beam_width);
  paras.Set<unsigned>("nodes_to_cache", 0);
  paras.Set<unsigned>("num_search_threads",
                      256);  // more for less contention of thread data.
  std::vector<T> data_load;
  size_t         dim{}, aligned_dim{};

  diskann::Timer timer;

  std::cout << "Loading queries " << std::endl;
  T*     query = NULL;
  size_t query_num, query_dim, query_aligned_dim;
  diskann::load_aligned_bin<T>(query_file, query, query_num, query_dim,
                               query_aligned_dim);

  dim = query_dim;
  aligned_dim = query_aligned_dim;
  diskann::Metric metric = diskann::Metric::L2;

  std::cout << "Searching before inserts: " << std::endl;

  int               batch = step * 10;  // step / 10;
  int               inMemorySize = 0;
  std::future<void> merge_future;

  fprintf(stderr, "Need copy: %d, need swap: %d\n", ckpt == 0, ckpt % 2 != 0);
  diskann::MergeInsert<T, TagT> sync_index(
      paras, dim, save_path + "_mem", save_path, save_path + "_merge", dist_cmp,
      metric, false, save_path, ckpt == 0, ckpt % 2 != 0, 0);

  uint64_t npts_modulo = sync_index._disk_index->num_points;
  uint64_t x = 1;
  while (npts_modulo > 1) {
    x *= 10;
    npts_modulo /= 10;
  }

  uint64_t index_npts = x;
  uint64_t vecs_per_step = index_npts / step;
  double   merge_ratio = ((double) Merge_Size / index_npts);

  uint64_t ckpt_i = 0, res = 0, tmp_npts = index_npts;
  for (int i = 0; i < ckpt; ++i) {
    res += Merge_Size;
    ckpt_i += Merge_Size / vecs_per_step;
    tmp_npts += Merge_Size;
    // Merge_Size =
    //     ((uint32_t) (merge_ratio * tmp_npts)) / vecs_per_step *
    //     vecs_per_step;
  }
  sync_index.init_mem_index(Merge_Size);
  LOG(INFO) << "index npts: " << index_npts
            << " vecs per step: " << vecs_per_step << " ckpt_i: " << ckpt_i
            << " merge ratio: " << merge_ratio << " res: " << res;
  LOG(INFO) << "Merge size is set to " << Merge_Size;
  // exit(0);

  if (ckpt_i == 0 && ckpt != 0) {
    std::cerr << "ckpt_i is 0, please check ckpt and Merge_Size" << std::endl;
    exit(-1);
  }

  std::string currentFileName = GetTruthFileName(truthset_file, res);
  begin_time = globalTimer.elapsed() / 1.0e6f;
  ShowMemoryStatus();
  for (int i = 0; i < Lsearch.size(); ++i) {
    sync_search_kernel(query, query_num, query_aligned_dim, recall_at,
                       Lsearch[i], sync_index, currentFileName, false, true);
  }

  for (int i = ckpt_i; i < batch; i++) {
    std::cout << "Batch: " << i << " Total Batch : " << batch << std::endl;
    std::vector<unsigned> insert_vec;
    std::vector<unsigned> delete_vec;

    /**Prepare for update*/
    uint64_t st = vecs_per_step * i;
    get_trace<T, TagT>(data_bin, st, st + index_npts, vecs_per_step, delete_vec,
                       insert_vec, data_load);

    std::future<void> insert_future = std::async(
        std::launch::async, insertion_kernel<T, TagT>, data_load.data(),
        std::ref(sync_index), std::ref(insert_vec), aligned_dim);

    std::future<void> delete_future = std::async(
        std::launch::async, deletion_kernel<T, TagT>, data_load.data(),
        std::ref(sync_index), std::ref(delete_vec), aligned_dim);

    int                total_queries = 0;
    std::future_status insert_status;
    std::future_status delete_status;
    do {
      insert_status = insert_future.wait_for(std::chrono::seconds(5));
      delete_status = delete_future.wait_for(std::chrono::seconds(5));
      if (insert_status == std::future_status::deferred ||
          delete_status == std::future_status::deferred) {
        std::cout << "deferred\n";
      } else if (insert_status == std::future_status::timeout ||
                 delete_status == std::future_status::timeout) {
        ShowMemoryStatus();
        sync_search_kernel(query, query_num, query_aligned_dim, recall_at,
                           Lsearch[0], sync_index, currentFileName, false,
                           false);
        total_queries += query_num;
        std::cout << "Queries processed: " << total_queries << std::endl;
      }
      if (insert_status == std::future_status::ready) {
        std::cout << "Insertions complete!\n";
      }
      if (delete_status == std::future_status::ready) {
        std::cout << "Deletions complete!\n";
      }
    } while (insert_status != std::future_status::ready ||
             delete_status != std::future_status::ready);

    inMemorySize += insert_vec.size();

    std::cout << "Search after update, current vector number: " << res
              << std::endl;

    res += vecs_per_step;
    currentFileName = GetTruthFileName(truthset_file, res);
    ShowMemoryStatus();

    for (int i = 0; i < Lsearch.size(); ++i) {
      sync_search_kernel(query, query_num, query_aligned_dim, recall_at,
                         Lsearch[i], sync_index, currentFileName, false, true);
    }

    if (i == batch - 1) {
      std::cout << "Done" << std::endl;
      exit(0);
    } else if (inMemorySize >= Merge_Size) {
      std::cout << "Begin Merge" << std::endl;
      merge_future = std::async(std::launch::async, merge_kernel<T, TagT>,
                                std::ref(sync_index), std::ref(save_path));
      std::this_thread::sleep_for(std::chrono::seconds(5));
      std::cout << "Sending Merge" << std::endl;
      inMemorySize = 0;
      std::future_status merge_status;
      do {
        merge_status = merge_future.wait_for(std::chrono::seconds(10));
        ShowMemoryStatus();
        sync_search_kernel(query, query_num, query_aligned_dim, recall_at,
                           Lsearch[0], sync_index, currentFileName, false,
                           false);
      } while (merge_status != std::future_status::ready);
      ShowMemoryStatus();
      std::cout << "Merge finished for checkpoint " << ckpt;
      exit(0);  // wait for reboot.
    }
  }
}

int main(int argc, char** argv) {
  if (argc < 15) {
    std::cout
        << "Correct usage: " << argv[0]
        << " <type[int8/uint8/float]> <data_bin> <L_disk> <R_disk> <alpha_disk>"
        << " <num_start> <#nodes_to_cache>"
        << " <indice_path> <query_file> <truthset_prefix> <recall@>"
        << " <#beam_width>"
        << " <step> <ckpt> <Lsearch> <L2>" << std::endl;
    exit(-1);
  }

  int         arg_no = 2;
  std::string data_bin = std::string(argv[arg_no++]);
  unsigned    L_disk = (unsigned) atoi(argv[arg_no++]);
  unsigned    R_disk = (unsigned) atoi(argv[arg_no++]);
  float       alpha_disk = (float) std::atof(argv[arg_no++]);
  size_t      num_start = (size_t) std::atoi(argv[arg_no++]);
  unsigned    nodes_to_cache = (unsigned) std::atoi(argv[arg_no++]);
  std::string save_path(argv[arg_no++]);

  std::string query_file(argv[arg_no++]);
  std::string truthset(argv[arg_no++]);
  int         recall_at = (int) std::atoi(argv[arg_no++]);
  unsigned    beam_width = (unsigned) std::atoi(argv[arg_no++]);
  int         step = (int) std::atoi(argv[arg_no++]);
  int         ckpt = (int) std::atoi(argv[arg_no++]);
  Merge_Size = (unsigned) std::atoi(argv[arg_no++]);
  std::cerr << "Merge size: " << Merge_Size << std::endl;

  std::vector<uint64_t> Lsearch;
  for (int i = arg_no; i < argc; ++i) {
    Lsearch.push_back(std::atoi(argv[i]));
  }
  for (auto& x : Lsearch) {
    std::cerr << "Lsearch: " << x << std::endl;
  }

  if (std::string(argv[1]) == std::string("int8")) {
    diskann::DistanceL2Int8 dist_cmp;
    update<int8_t, unsigned>(data_bin, L_disk, R_disk, alpha_disk, step,
                             num_start, nodes_to_cache, save_path, query_file,
                             truthset, recall_at, Lsearch, beam_width,
                             &dist_cmp, ckpt);
  } else if (std::string(argv[1]) == std::string("uint8")) {
    diskann::DistanceL2UInt8 dist_cmp;
    update<uint8_t, unsigned>(data_bin, L_disk, R_disk, alpha_disk, step,
                              num_start, nodes_to_cache, save_path, query_file,
                              truthset, recall_at, Lsearch, beam_width,
                              &dist_cmp, ckpt);
  } else if (std::string(argv[1]) == std::string("float")) {
    diskann::DistanceL2 dist_cmp;
    update<float, unsigned>(data_bin, L_disk, R_disk, alpha_disk, step,
                            num_start, nodes_to_cache, save_path, query_file,
                            truthset, recall_at, Lsearch, beam_width, &dist_cmp,
                            ckpt);
  } else
    std::cout << "Unsupported type. Use float/int8/uint8" << std::endl;
}
