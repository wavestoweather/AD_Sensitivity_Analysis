#include "include/types/global_args_t.h"

global_args_t::global_args_t()
{
    final_time_flag = 0;
    final_time_string = nullptr;

    timestep_flag = 0;
    timestep_string = nullptr;

    snapshot_index_flag = 0;
    snapshot_index_string = nullptr;

    output_flag = 0;
    output_string = nullptr;

    scaling_fact_flag = 0;
    scaling_fact_string = nullptr;

    input_flag = 0;
    input_file = nullptr;

    start_over_flag = 0;
    start_over_string = nullptr;

    start_over_env_flag = 0;
    start_over_env_string = nullptr;

    fixed_iteration_flag = 0;
    fixed_iteration_string = nullptr;

    auto_type_flag = 0;
    auto_type_string = nullptr;

    traj_flag = 0;
    traj_string = nullptr;

    write_flag = 0;
    write_string = nullptr;

    progress_index_flag = 0;
    progress_index_string = nullptr;

#ifdef MET3D
    delay_start_flag = 0;
    delay_start_string = nullptr;
#endif

    ens_config_flag = 0;
    ens_config_string = nullptr;

    checkpoint_flag = 0;
    checkpoint_string = nullptr;

    gnu_id_flag = 0;
    gnu_id_string = nullptr;

    folder_name_flag = 0;
    folder_name_string = nullptr;

    n_ens_flag = 0;
    n_ens_string = nullptr;
}

int global_args_t::parse_arguments(
    const int argc,
    char* const * argv,
    const int &rank,
    const int &n_processes)
{
    /**
     * String used to parse commandline input.
     */
    static const char *optString = "w:f:d:e:i:b:o:l:s:t:a:r:p:n:m:c:g:h:k:?";
    bool need_to_abort = false;
    int opt;

    if(argc < 2)
    {
        need_to_abort = true;
        display_usage();
    }else
    {
        opt = getopt(argc, argv, optString);

        while(-1 != opt)
        {
            switch(opt)
            {
                case 'f':
                {
                    this->final_time_flag = 1;
                    this->final_time_string = optarg;
                    break;
                }
                case 'd':
                {
                    this->timestep_flag = 1;
                    this->timestep_string = optarg;
                    break;
                }
                case 'i':
                {
                    this->snapshot_index_flag = 1;
                    this->snapshot_index_string = optarg;
                    break;
                }
                case 'b':
                {
                    this->scaling_fact_flag = 1;
                    this->scaling_fact_string = optarg;
                    break;
                }
                case 'o':
                {
                    this->output_flag = 1;
                    this->output_string = optarg;
                    break;
                }
                case 'l':
                {
                    this->input_flag = 1;
                    this->input_file = optarg;
                    break;
                }
                case 's':
                {
                    this->start_over_flag = 1;
                    this->start_over_string = optarg;
                    break;
                }
                case 'e':
                {
                    this->start_over_env_flag = 1;
                    this->start_over_env_string = optarg;
                    break;
                }
                case 't':
                {
                    this->fixed_iteration_flag = 1;
                    this->fixed_iteration_string = optarg;
                    break;
                }
                case 'a':
                {
                    this->auto_type_flag = 1;
                    this->auto_type_string = optarg;
                    break;
                }
                case 'r':
                {
                    this->traj_flag = 1;
                    this->traj_string = optarg;
                    break;
                }
                case 'w':
                {
                    this->write_flag = 1;
                    this->write_string = optarg;
                    break;
                }
                case 'p':
                {
                    this->progress_index_flag = 1;
                    this->progress_index_string = optarg;
                    break;
                }
#ifdef MET3D
                case 'n':
                {
                    this->delay_start_flag = 1;
                    this->delay_start_string = optarg;
                    break;
                }
#endif
                case 'm':
                {
                    this->ens_config_flag = 1;
                    this->ens_config_string = optarg;
                    break;
                }
                case 'c':
                {
                    this->checkpoint_flag = 1;
                    this->checkpoint_string = optarg;
                    break;
                }
                case 'g':
                {
                    this->gnu_id_flag = 1;
                    this->gnu_id_string = optarg;
                    break;
                }
                case 'h':
                {
                    this->folder_name_flag = 1;
                    this->folder_name_string = optarg;
                    break;
                }
                case 'k':
                {
                    this->n_ens_flag = 1;
                    this->n_ens_string = optarg;
                    break;
                }
                case '?':
                {
                    need_to_abort = true;
                    display_usage();
                    break;
                }
                default:
                {
                    need_to_abort = true;
                    display_error_on_command_line();
                    display_usage();
                    break;
                }
            }

            opt = getopt(argc, argv, optString);
        }
    }

    if(need_to_abort){
        return ARGUMENT_ERR;
    }
    return SUCCESS;
}

void global_args_t::display_usage()
{
    std::cout << "\n"
        << "USAGE of the program:\n"
        << "Invoke the program on the command line with\n"
        << "$ ./main PARAMETERS\n"
        << "where PARAMETERS are the following parameters:\n"
        << "-f: Time to integrate in seconds. "
        << "If higher than number of entries from input file, then it stops early\n"
        << "-i: Snapshot index.\n"
        << "-b: Scaling factor.\n"
        << "-d: Timestep in seconds for a substep between each new trajectory input.\n"
        << "-o: Name of the output file.\n"
        << "-l: Path and name of input file.\n"
        << "-s: Start over mixing ratios and particle numbers at each timestep of a trajectory.\n"
        << "-e: Start over temperature, pressure and ascent at each timestep of a trajectory.\n"
        << "-t: Set pressure, temperature and vertical velocity fixed at each substep.\n"
        << "-a: Set auto type (1=KB, 2=KK, 3=SB).\n"
        << "-r: Set index of trajectory.\n"
        << "-w: Write index for the snapshots.\n"
        << "-p: Index for updating the progressbar.\n"
#ifdef MET3D
        << "-n: No simulation until given time (relative to ascend).\n"
#endif
        << "-m: Name of an ensemble configuration file.\n"
        << "-c: Name of a checkpoint file.\n"
        << "-g: ID (uint) for this instance.\n"
        << "-h: Path and name to subsequent folders for any checkpoints and "
        << "execution scripts from this simulation.\n"
        << "-k: Maximum number of ensembles in the output file. If none given "
        << "then each segment from the configuration file is considered as one ensemble.\n"
        << "-?: This help message.\n"
        << std::endl;
}

void global_args_t::display_error_on_command_line()
{
    std::cerr << "==> ERROR: An error occured while dealing with the command line arguments!"
        << std::endl;
}