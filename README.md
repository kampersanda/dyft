# Dynamic Filter Trie (DyFT)

This is an experimental library for data structures described in the following papers:

- Shunsuke Kanda and Yasuo Tabei, [Dynamic Similarity Search on Integer Sketches](https://ieeexplore.ieee.org/document/9338383/), In *Proc. 20th IEEE ICDM*, pp 242–251, 2020 [[arXiv](https://arxiv.org/abs/2009.11559)]
- Shunsuke Kanda and Yasuo Tabei, [DyFT: A Dynamic Similarity Search Method on Integer Sketches](https://link.springer.com/article/10.1007%2Fs10115-021-01611-2), *Knowledge and Information Systems*, 63, 2815–2840, 2021 [SharedIt](https://rdcu.be/cxu1J)

## Build instructions

You can download and compile this library by the following commands:

```shell
$ git clone https://github.com/kampersanda/dyft.git
$ cd dyft
$ mkdir build
$ cd build
$ cmake ..
$ make -j
```

After the commands, the executables will be produced in `build/bin` directory.

The code is written in C++17, so please install g++ >= 7.0 or clang >= 4.0. The following dependencies have to be installed to compile the library: CMake >= 3.0, for the build system, and Boost >= 1.42.

The library employs the third-party libraries [cmd\_line\_parser](https://github.com/jermp/cmd_line_parser) and [tinyformat](https://github.com/c42f/tinyformat), whose header files are contained in this repository.

## Data structures

The library contains some implementations of data structures in addition to our DyFT, which were used in the experiments (Section VI of our paper), as follows.

### DyFT family

- `mart_index` is an implementation of DyFT with MART representation.
- `art_index` is an implementation of DyFT with ART representation.
- `array_index` is an implementation of DyFT with Array representation.
- `mi_frame<T>` is an implementation of DyFT+ (i.e., multi-index variant) for the index which is specified by the template argument `T` of a DyFT class.

### Others

- `hms1v_index` is an implementation of [HmSearch 1-var (HSV)](https://doi.org/10.1145/2484838.2484842), which is designed for binary sketches.
- `hms1dv_index` is an implementation of [HmSearch 1-del-var (HSD)](https://doi.org/10.1145/2484838.2484842), which is designed for integer sketches.
- `gv_index` is an implementation of multi-index hashing for integer sketches, which follows an idea by [Gog and Venturini](https://doi.org/10.1145/2911451.2911523).

## Input data format

### Binary sketches

Binary sketches should be stored in binary format, where each sketch is 64 bits of size. The executables provided by the library can indicate the number of dimensions (32 or 64) to evaluate through argument `-N`, and use the first `N` dimensions.

You can use `src/simhash.cpp` to generate such a dataset using [Charikar’s simhash](https://doi.org/10.1145/509907.509965) algorithm.

### Integer sketches

Integer sketches should be stored in [TEXMEX's bvecs format](http://corpus-texmex.irisa.fr/). That is, the dimension number and features for each sketch are interleaved, where the number is 4 bytes of size and each feature is 1 byte of size. In the same manner, the executables provided by the library can indicate the dimension number (32 or 64) to evaluate through argument `-N`. And, argument `-B` can indicate the number of lowest bits to evaluate for each feature (1 to 8).

You can use [consistent\_weighted\_sampling](https://github.com/kampersanda/consistent_weighted_sampling) to generate such a dataset using the [GCWS](https://doi.org/10.1145/3097983.3098081) algorithm.

## Usage

In the `data` directory, there are the tiny datasets:

- `SIFT.10K.bin` and `SIFT.100.bin` are binary sketches generated from [ANN\_SIFT10K](http://corpus-texmex.irisa.fr/) using Charikar’s simhash algorithm.
- `mnist.scale.cws.bvecs` and `mnist.scale.1K.cws.bvecs` are 1-byte integer sketches generated from [minst](https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multiclass.html#mnist) using the GCWS algorithm.

By using these datasets, you can test the data structure as follows.

### Benchmark for range search

The executables `range_search_*` are used to analyze the performance of range search time for each data structure. Given a dataset and query files, the executable inserts data sketches one by one and performs range search for the queries when the number of inserted sketches is power of ten. The benchmark results are output as a json file.

#### Example: DyFT with MART on binary sketches

```shell
$ ./bin/range_search_bin_dyft ../data/SIFT.10K.bin ../data/SIFT.100.bin -o results -R 2 -N 32 -A mart -K 1
```

The command tests range search with radius `R=2` for `mart_index` of single index (`K=1`) through dataset `SIFT.10K.bin` and queryset `SIFT.100.bin`, where the number of dimensions is `N=32` (i.e., the first 32 dimensions are used). The benchmark result is output in `results` directory, where the file name is assigned automatically.

#### Example: DyFT+ with ART on integer sketches

```shell
$ ./bin/range_search_int_dyft ../data/mnist.scale.cws.bvecs ../data/mnist.scale.1K.cws.bvecs -o results -R 4 -N 64 -B 4 -A art -K 3
```

The command tests range search with radius `R=4` for `mi_frame<art_index>` of `K=3` blocks through dataset `mnist.scale.cws.bvecs` and queryset `mnist.scale.1K.cws.bvecs`, where the number of dimensions is `N=64` and the lowest `B=4` bits are used for each feature.

#### Example: GV on integer sketches

```shell
$ ./bin/range_search_int_gv ../data/mnist.scale.cws.bvecs ../data/mnist.scale.1K.cws.bvecs -o results -R 6 -N 64 -B 8
```

The command tests range search with radius `R=6` for `gv_index` through dataset `mnist.scale.cws.bvecs` and queryset `mnist.scale.1K.cws.bvecs`, where the number of dimensions is `N=64` and all the `B=8` bits are used for each feature.

### Benchmark for construction

The executables `build_index_*` are used to analyze the performance of insertion time and process size for each data structure. Given a dataset, the executable inserts data sketches one by one and reports the statistics when the number of inserted sketches is power of ten. The benchmark results are output as a json file.

#### Example: DyFT with Array on binary sketches

```shell
$ ./bin/build_index_bin_dyft ../data/SIFT.10K.bin -o results -R 2 -N 32 -A array -K 1
```

The command builds `array_index` of single index (`K=1`) on radius `R=2` for dataset `SIFT.10K.bin`, where the number of dimension is `N=32`.

#### Example: DyFT+ with MART on integer sketches

```shell
$ ./bin/build_index_int_dyft ../data/mnist.scale.cws.bvecs -o results -R 4 -N 64 -B 4 -A mart -K 0
```

The command builds `mi_frame<mart_index>` with a reasonable number of blocks on radius `R=4` for dataset `mnist.scale.cws.bvecs`, where the number of dimension is `N=64` and the lowest `B=4` bits are used for each feature.

When `K` is set to 0, the number of blocks is set to ⌊`R`/2⌋+1 following the idea of GV.

#### Example: HSV on binary sketches

```shell
$ ./bin/build_index_bin_hms1v ../data/SIFT.10K.bin -o results -R 6 -N 64
```

The command builds `hms1v_index` on radius `R=6` for dataset `SIFT.10K.bin`, where the number of dimension is `N=64`.

## Licensing

This library is free software provided under [MIT License](https://github.com/kampersanda/dyft/blob/master/LICENSE).

If you use the library, please cite the following paper:

```
@inproceedings{kanda2020dynamic,
  author = {Kanda, Shunsuke and Tabei, Yasuo},
  title = {Dynamic Similarity Search on Integer Sketches},
  booktitle = {Proceedings of the 20th IEEE International Conference on Data Mining (ICDM)},
  pages={242-251},
  year = {2020}
}
```

