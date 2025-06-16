# loop from 1 to 100
for i in {0..18}
do
    # Path: build/overall.sh
    # Run the build script
    echo "Running with i = $i"
    OMP_PLACES=cores OMP_PROC_BIND=close build/tests/motivation int8 /mnt/nvme/data/SPACEV1B/vectors.bin 128 64 1.2 0 0 /mnt/nvme/indices/SPACEV1B/100M /mnt/nvme/data/SPACEV1B/query.bin /mnt/nvme/data/SPACEV1B/gnd_insert/200M_topk 10 4 100 $i 6000000 20 25 30 33 37 40
    # build/tests/motivation uint8 /mnt/nvme/data/bigann/bigann_2M.bbin 128 64 1.2 0 0 /mnt/nvme/indices/bigann/1m /mnt/nvme/data/bigann/bigann_query.bbin /mnt/nvme/data/bigann/gnd/2M_topk 10 4 100 $i 20000 20 30 40
done
