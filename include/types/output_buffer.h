#pragma once

#include <array>
#include <string>
#include <vector>

struct output_buffer_vector_t{
    // for netCDF files and a vector for each column
    // columns: output parameter + gradients + time_afer_ascent + type + flags
    // fast index: record
    // each array = one column
    // slow index: num_comp
    std::array<std::vector<double>, num_comp+num_par+4 > output_buffer;
    std::array<std::vector<unsigned char>, 4 > output_buffer_flags;
    std::array<std::vector<std::string>, 1 > output_buffer_str;
    std::array<std::vector<uint64_t>, 1 > output_buffer_int;

    std::array<double,
    uint64_t n_snapshots; // number of buffered snapshots
}

// Idea: Every new checkpoint gets traj_id = rank + i*n_processes where i is
// the number of processes checkpoints by the current process
struct output_max_buffer_t{
    // for netCDF files and a vector for each column
    // columns: output parameter + gradients + time_afer_ascent + type + flags
    // fast index: record
    // each array = one column
    // slow index: num_comp
    uint64_t n_snapshots;
    std::array<double, 1400*(num_comp + num_par*num_comp + 4 ) > output_buffer;
    std::array<unsigned char, 1400*4 > output_buffer_flags;
    std::array<std::vector<std::string>, 1400 > output_buffer_str;
    std::array<std::vector<uint64_t>, 1400 > output_buffer_int;
}

struct output_buffer_t{
    // for netCDF files and a vector for each column
    // columns: output parameter + gradients + time_afer_ascent + type + flags
    // fast index: record
    // each array = one column
    // slow index: num_comp
    std::array<double, (num_comp + num_par*num_comp + 4 ) > output_buffer;
    std::array<unsigned char, 4 > output_buffer_flags;
    std::array<std::vector<std::string>, 1 > output_buffer_str;
    std::array<std::vector<uint64_t>, 1 > output_buffer_int;
}

MPI_Datatype MPI_OUTPUT_BUFFER;
MPI_Datatype MPI_MAX_OUTPUT_BUFFER;

void create_MPI_type();