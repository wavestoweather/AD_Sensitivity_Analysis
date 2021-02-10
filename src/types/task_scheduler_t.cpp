#include "include/types/task_scheduler_t.h"

task_scheduler_t::task_scheduler_t(
    const int &rank,
    const int &n_processes)
{
    free_worker.resize(n_processes);
    if(rank == 0)
    {
        std::fill(free_worker.begin()+1, free_worker.end(), 1);
        free_worker[0] = 0;
    } else
    {
        std::fill(free_worker.begin(), free_worker.end(), 0);
    }

    work_available.resize(n_processes);
    std::fill(work_available.begin(), work_available.end(), 0);
    SUCCESS_OR_DIE(
        MPI_Win_create(
            free_worker.data(),
            n_processes*sizeof(std::int8_t),
            sizeof(std::int8_t),
            MPI_INFO_NULL,
            MPI_COMM_WORLD,
            &free_window)
    );
    SUCCESS_OR_DIE(
        MPI_Win_create(
            work_available.data(),
            n_processes*sizeof(std::int8_t),
            sizeof(std::int8_t),
            MPI_INFO_NULL,
            MPI_COMM_WORLD,
            &work_window)
    );

    my_rank = rank;
}

void task_scheduler_t::send_task()
{
    // Send a task from own queue if some worker is free
    if(!checkpoint_queue.empty())
    {
        auto it = std::find(free_worker.begin(), free_worker.end(), 1);
        if(it != free_worker.end())
        { // Send task from own queue to someone
            checkpoint_t checkpoint = checkpoint_queue.front();
            checkpoint_queue.pop();
            checkpoint.send_checkpoint(*it);
            if(my_rank == 0)
            {
                free_worker[*it] = 0;
                work_available[0]--;
            }
        }
    }
    signal_send_task();
}

void task_scheduler_t::send_task(
    checkpoint_t &checkpoint)
{
    auto it = std::find(free_worker.begin(), free_worker.end(), 1);
    if(it != free_worker.end())
    {
        checkpoint.send_checkpoint(*it);
        free_worker[*it] = 0;
    } else
    {
        checkpoint_queue.push(checkpoint);
        work_available[my_rank]++;
        if(my_rank != 0)
        {
            signal_work_avail();
        }
    }

}

void task_scheduler_t::signal_send_task()
{
    // Check for any remaining tasks on any other process and signal
    // them where to send it.
    if(my_rank == 0)
    {
        uint32_t work_idx = 0;
        for(auto &worker: free_worker)
        {
            if(worker == 1)
            {
                for(uint32_t i=work_idx; i<work_available.size(); ++i)
                {
                    if(work_available[i] > 0)
                    {
                        // signal to work_idx to send to worker
                        std::int8_t free = 1;
                        SUCCESS_OR_DIE(
                            MPI_Put(
                                &free, // source
                                1, // count
                                MPI_INT8_T, // datatype
                                work_idx, // target
                                worker, // displacement in target
                                1, // target count
                                MPI_INT8_T, // target datatype
                                free_window)
                        );

                        // Mark worker as busy and remove work
                        free_worker[worker] = 0;
                        work_available[i]--;
                        work_idx = i;
                        break; // continue with next free worker
                    }
                }
            }
        }
    }
}

bool task_scheduler_t::receive_task(
    checkpoint_t &checkpoint)
{
    bool signal_sent = false;
    // check if work is directly available
    while(!checkpoint.receive_checkpoint())
    {
        // if not, signal free to the master and continue busy waiting
        if(!signal_sent && my_rank != 0)
        {
            signal_free();
            signal_sent = true;
        }

        // if all are free, return
        if(my_rank == 0)
        {
            auto it = std::find(free_worker.begin(), free_worker.end(), 1);
            bool no_more_work = true;
            for(const auto &work: work_available)
            {
                if(work > 0)
                {
                    no_more_work = false;
                    break;
                }
            }
            if(it == free_worker.end() && no_more_work)
            {
                // signal everyone that we are finished somehow
                free_worker[0] = 2;
                // We send it to every worker, however broadcast a single
                // value everytime we are here is not a good idea.
                for(uint32_t i=0; i<free_worker.size(); ++i)
                {
                    SUCCESS_OR_DIE(
                    MPI_Put(
                        free_worker.data(), // source
                        1, // count
                        MPI_INT8_T, // datatype
                        i, // target
                        0, // displacement in target
                        1, // target count
                        MPI_INT8_T, // target datatype
                        free_window)
                    );
                }
                return false;
            }
        } else if(free_worker[0] == 2)
        {
            return false;
        }

    };
    return true;
}

void task_scheduler_t::signal_free()
{
    std::int8_t free = 1;
    SUCCESS_OR_DIE(
        MPI_Put(
            &free, // source
            1, // count
            MPI_INT8_T, // datatype
            0, // target
            my_rank, // displacement in target
            1, // target count
            MPI_INT8_T, // target datatype
            free_window)
    );
    // MPI_Request request;
    // free_worker[my_rank] = 1;
    // SUCCESS_OR_DIE(
    //     MPI_Isend(
    //         &free_worker[my_rank],
    //         1,
    //         MPI_INT,
    //         my_rank,
    //         MPI_COMM_WORLD,
    //         &request)
    // );
}

void task_scheduler_t::signal_work_avail()
{
    std::int8_t work = checkpoint_queue.size();
    SUCCESS_OR_DIE(
        MPI_Put(
            &work, // source
            1, // count
            MPI_INT8_T, // datatype
            0, // target
            my_rank, // displacement in target
            1, // target count
            MPI_INT8_T, // target datatype
            work_window)
    );
}