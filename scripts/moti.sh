# loop from 1 to 100
for i in 1000000 2000000 4000000 6000000 8000000 16000000 32000000
do
    # Path: build/overall.sh
    # Run the build script
    echo "Running with i = $i"
    OMP_PROC_BIND=close build/tests/motivation uint8 /mnt/nvme/data/bigann/bigann.bbin 128 64 1.2 0 0 /mnt/nvme/indices/bigann/100m /mnt/nvme/data/bigann/bigann_query.bbin /mnt/nvme/data/bigann/gnd/500M_topk 10 4 100 0 $i 20 30 40
    # build/tests/motivation uint8 /mnt/nvme/data/bigann/bigann_2M.bbin 128 64 1.2 0 0 /mnt/nvme/indices/bigann/1m /mnt/nvme/data/bigann/bigann_query.bbin /mnt/nvme/data/bigann/gnd/2M_topk 10 4 100 0 10000 20 30 40
done
