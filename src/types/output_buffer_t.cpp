#include "include/types/output_buffer_t.h"

void create_MPI_type()
{
    std::vector<MPI_Datatype> data_types(4) = {
        MPI_DOUBLE, MPI_CHAR, MPI_CHAR, MPI_UINT64_T
    };
    std::vector<int> len(4);
    std::vector<int> disp(4);
    len[0] = (num_comp + num_par*num_comp + 4 );
    len[1] = 4;
    len[2] = 1;
    len[3] = 1;

    disp[0] = 0;
    disp[1] = disp[0] + len[0]*std::sizeof(double);
    disp[2] = disp[1] + len[1]*std::sizeof(char);
    disp[3] = disp[2] + len[2]*std::sizeof(char);

    MPI_Datatype tmp;
    SUCCESS_OR_DIE(
        MPI_Type_create_struct(
        5, len.data(), disp.data(), data_types.data(), &tmp)
    );
    SUCCESS_OR_DIE(
        MPI_TYPE_create_resized(tmp, 0, sizeof(output_buffer_t),
        &MPI_OUTPUT_BUFFER)
    );
    SUCCESS_OR_DIE(
        MPI_TYPE_commit(&MPI_OUTPUT_BUFFER)
    );
    SUCCESS_OR_DIE(
        MPI_Type_free(&tmp)
    );
}
void create_MPI_type_max()
{
    std::vector<MPI_Datatype> data_types(5) = {
        MPI_UINT64_T, MPI_DOUBLE, MPI_CHAR, MPI_CHAR, MPI_UINT64_T
    };
    std::vector<int> len(5);
    std::vector<int> disp(5);
    // Same with maxed out buffer
    len[0] = 1;
    len[1] = 1400*(num_comp + num_par*num_comp + 4 );
    len[2] = 1400*4;
    len[3] = 1400;
    len[4] = 1400;

    disp[0] = 0;
    disp[1] = std::sizeof(uint64_t);
    disp[2] = disp[1] + len[1]*std::sizeof(double);
    disp[3] = disp[2] + len[2]*std::sizeof(char);
    disp[4] = disp[3] + len[3]*std::sizeof(char);

    MPI_Datatype tmp;
    SUCCESS_OR_DIE(
        MPI_Type_create_struct(
        5, len.data(), disp.data(), data_types.data(), &tmp)
    );
    SUCCESS_OR_DIE(
        MPI_TYPE_create_resized(tmp, 0, sizeof(output_max_buffer_t),
        &MPI_MAX_OUTPUT_BUFFER)
    );
    SUCCESS_OR_DIE(
        MPI_TYPE_commit(&MPI_MAX_OUTPUT_BUFFER)
    );
    SUCCESS_OR_DIE(
        MPI_Type_free(&tmp)
    );
}